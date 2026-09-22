#!/usr/bin/env python3
"""
annotate_variants.py

Annotate breast-cancer risk variants (RUNFILE + ExpectoSC/Seqweaver z-score
output) with genomic context - overlapping gene(s), exonic/intronic/UTR
location, and position relative to the TSS - using a GTF gene annotation,
then attach the per-cell-type z-scores and a few summary statistics.

------------------------------------------------------------------------
GENOME BUILD
------------------------------------------------------------------------
The RUNFILE coordinates are GRCh37/hg19 (verified: the RUNFILE places the
HPDL TSS at chr1:45,792,568, which matches HPDL's GRCh37.p13 RefSeq
annotation at chr1:45,792,567-45,794,382; on GRCh38 that gene sits at
chr1:45,326,905 - a completely different coordinate, ruling that build
out). Use an hg19/GRCh37 GTF. Good sources:

  Ensembl GRCh37 (recommended - bare chromosome names "1","2",... like
  your RUNFILE, so no chr-prefix stripping needed, though this script
  handles either convention):
    https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

  GENCODE lift37 (uses "chr1" naming):
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh37_mapping/gencode.v48lift37.annotation.gtf.gz

------------------------------------------------------------------------
INPUT FILE ASSUMPTIONS
------------------------------------------------------------------------
RUNFILE columns (as seen in Seqweaver_bridges34_genes_input_RUNFILE.csv):
    chrom, snp_chromStart, snp_chromEnd, SNP, ref, alt, chrom_gene,
    tss, tss_chromEnd, strand, gene_name, gene_id, cage_used, dist
  - one row per (variant, tested-gene) pair - a variant can appear more
    than once if it was screened against more than one gene
  - snp_chromStart/snp_chromEnd look 0-based/half-open (UCSC-BED style);
    this script uses snp_chromEnd as the 1-based genomic position,
    matching the tss/tss_chromEnd convention in the same file.

Z-score file columns (Seqweaver_bridges__breast.csv):
    SNP, <cell_type_1>, <cell_type_2>, ... <cell_type_59>
  - one row per unique variant, one column per cell type

------------------------------------------------------------------------
OUTPUT
------------------------------------------------------------------------
One row per (variant, tested-gene) pair with:
  - variant/gene/TSS metadata from the RUNFILE
  - genes_overlapping_variant : any gene(s) from the GTF whose body
    (gene feature start-end) contains the variant - independent of
    which gene it was tested against, useful as a QC check
  - feature_type   : exonic_CDS / exonic_UTR / exonic / intronic /
                      intergenic, relative to the TESTED gene
                      (gene_id column) specifically
  - overlapping_transcripts : transcript IDs of the tested gene whose
                      exon the variant falls in (blank if none)
  - position_vs_tss : upstream / downstream / at_tss, strand-aware
  - max_abs_z, max_abs_z_celltype, mean_z, n_celltypes_abs_z_ge_<thr>
  - all per-cell-type z-score columns

------------------------------------------------------------------------
USAGE
------------------------------------------------------------------------
    python3 annotate_variants.py \
        --runfile Seqweaver_bridges34_genes_input_RUNFILE.csv \
        --zscores Seqweaver_bridges__breast.csv \
        --gtf Homo_sapiens.GRCh37.87.gtf.gz \
        --out annotated_variants.csv \
        --z-threshold 2.0
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree

# ---------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------

# Feature types we care about. GENCODE (newer releases) uses
# five_prime_utr / three_prime_utr; Ensembl and older GENCODE use "UTR".
EXON_FEATURES = {"exon"}
CDS_FEATURES = {"CDS"}
UTR_FEATURES = {"UTR", "five_prime_utr", "three_prime_utr"}
GENE_FEATURES = {"gene"}
KEEP_FEATURES = EXON_FEATURES | CDS_FEATURES | UTR_FEATURES | GENE_FEATURES

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def _open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def strip_version(gene_id):
    """ENSG00000186603.7 -> ENSG00000186603"""
    return gene_id.split(".")[0] if gene_id else gene_id


def load_gtf_features(gtf_path, target_gene_ids, target_gene_names):
    """
    Single pass over the GTF, keeping only gene/exon/CDS/UTR rows that
    belong to one of the target genes (matched by versionless gene_id
    OR gene_name, since RUNFILE gene_id may or may not carry the GTF's
    version suffix).

    Returns:
        gene_span      : dict gene_key -> (chrom, start, end)
        exon_trees      : dict gene_key -> IntervalTree of (start,end,transcript_id)
        cds_trees       : dict gene_key -> IntervalTree
        utr_trees       : dict gene_key -> IntervalTree
        gene_key_lookup : dict (versionless gene_id) -> gene_key,
                           dict gene_name -> gene_key
    All coordinates kept 1-based inclusive, as GTF defines them.
    """
    target_gene_ids = {strip_version(g) for g in target_gene_ids}
    target_gene_names = set(target_gene_names)

    gene_span = {}
    exon_trees = defaultdict(IntervalTree)
    cds_trees = defaultdict(IntervalTree)
    utr_trees = defaultdict(IntervalTree)
    id_to_key = {}
    name_to_key = {}

    n_lines = 0
    n_kept = 0
    with _open_maybe_gz(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_lines += 1
            # cheap pre-filter before full split/regex
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in KEEP_FEATURES:
                continue

            attr = fields[8]
            gid_m = GENE_ID_RE.search(attr)
            if not gid_m:
                continue
            raw_gene_id = gid_m.group(1)
            gene_id_nov = strip_version(raw_gene_id)

            gname_m = GENE_NAME_RE.search(attr)
            gene_name = gname_m.group(1) if gname_m else None

            if gene_id_nov not in target_gene_ids and gene_name not in target_gene_names:
                continue

            n_kept += 1
            chrom = fields[0].replace("chr", "")
            start = int(fields[3])
            end = int(fields[4])
            gene_key = gene_id_nov  # canonical key

            id_to_key[gene_id_nov] = gene_key
            if gene_name:
                name_to_key[gene_name] = gene_key

            if feature in GENE_FEATURES:
                gene_span[gene_key] = (chrom, start, end)
            elif feature in EXON_FEATURES:
                tid_m = TRANSCRIPT_ID_RE.search(attr)
                tid = tid_m.group(1) if tid_m else ""
                if end >= start:
                    exon_trees[gene_key][start:end + 1] = tid
            elif feature in CDS_FEATURES:
                if end >= start:
                    cds_trees[gene_key][start:end + 1] = True
            elif feature in UTR_FEATURES:
                if end >= start:
                    utr_trees[gene_key][start:end + 1] = True

    print(
        f"[gtf] scanned {n_lines:,} lines, kept {n_kept:,} rows for "
        f"{len(gene_span)} target genes with a 'gene' record "
        f"(target list had {len(target_gene_ids)} gene_ids / "
        f"{len(target_gene_names)} gene_names)",
        file=sys.stderr,
    )
    missing = target_gene_ids - set(gene_span.keys())
    if missing:
        print(
            f"[gtf] WARNING: {len(missing)} target gene_id(s) not found as a "
            f"'gene' feature in the GTF (check chromosome naming / gene_id "
            f"version / GTF build): {sorted(missing)[:10]}"
            + (" ..." if len(missing) > 10 else ""),
            file=sys.stderr,
        )

    return gene_span, exon_trees, cds_trees, utr_trees, id_to_key, name_to_key


# ---------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------

def classify_position(gene_key, pos, gene_span, exon_trees, cds_trees, utr_trees):
    """
    Returns (feature_type, overlapping_transcript_ids_str)
    feature_type in: exonic_CDS, exonic_UTR, exonic, intronic, intergenic,
                      unknown_gene (gene not found in GTF at all)
    """
    if gene_key not in gene_span:
        return "unknown_gene", ""

    _, gstart, gend = gene_span[gene_key]
    if not (gstart <= pos <= gend):
        return "intergenic", ""

    exon_hits = exon_trees[gene_key][pos] if gene_key in exon_trees else set()
    if not exon_hits:
        return "intronic", ""

    transcript_ids = ";".join(sorted({iv.data for iv in exon_hits if iv.data}))

    if gene_key in cds_trees and cds_trees[gene_key][pos]:
        return "exonic_CDS", transcript_ids
    if gene_key in utr_trees and utr_trees[gene_key][pos]:
        return "exonic_UTR", transcript_ids
    return "exonic", transcript_ids


def genes_overlapping(chrom, pos, gene_span):
    """All target genes (any of them, not just the tested one) whose gene
    body spans this position - used as an independent QC cross-check."""
    hits = [
        gk for gk, (gchrom, gstart, gend) in gene_span.items()
        if gchrom == chrom and gstart <= pos <= gend
    ]
    return hits


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runfile", required=True, help="RUNFILE CSV (variant/gene/TSS metadata)")
    ap.add_argument("--zscores", required=True, help="ExpectoSC/Seqweaver z-score CSV (one row per variant)")
    ap.add_argument("--gtf", required=True, help="hg19/GRCh37 GTF file (.gtf or .gtf.gz)")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument(
        "--z-threshold", type=float, default=2.0,
        help="abs(z) threshold used for the 'n_celltypes_abs_z_ge_<thr>' summary column (default: 2.0)",
    )
    args = ap.parse_args()

    print("[load] reading RUNFILE...", file=sys.stderr)
    run = pd.read_csv(args.runfile, dtype={"chrom": str, "chrom_gene": str})

    print("[load] reading z-score file...", file=sys.stderr)
    z = pd.read_csv(args.zscores)
    celltype_cols = [c for c in z.columns if c != "SNP"]

    missing_z = set(run.SNP) - set(z.SNP)
    missing_run = set(z.SNP) - set(run.SNP)
    if missing_z:
        print(f"[warn] {len(missing_z)} RUNFILE variants have no z-score row (will be dropped)", file=sys.stderr)
    if missing_run:
        print(f"[warn] {len(missing_run)} z-score variants are not in the RUNFILE (ignored)", file=sys.stderr)

    # ---- z-score summary stats, computed once per unique variant ----
    print("[compute] summarizing z-scores per variant...", file=sys.stderr)
    zvals = z[celltype_cols]
    abs_zvals = zvals.abs()
    z_summary = pd.DataFrame({
        "SNP": z["SNP"],
        "max_abs_z": abs_zvals.max(axis=1),
        "max_abs_z_celltype": abs_zvals.idxmax(axis=1),
        "mean_z": zvals.mean(axis=1),
        f"n_celltypes_abs_z_ge_{args.z_threshold:g}": (abs_zvals >= args.z_threshold).sum(axis=1),
    })

    # ---- GTF ----
    target_gene_ids = run["gene_id"].dropna().unique().tolist()
    target_gene_names = run["gene_name"].dropna().unique().tolist()
    print(f"[gtf] loading {args.gtf} (filtering to {len(target_gene_ids)} target genes)...", file=sys.stderr)
    gene_span, exon_trees, cds_trees, utr_trees, id_to_key, name_to_key = load_gtf_features(
        args.gtf, target_gene_ids, target_gene_names
    )

    def resolve_gene_key(gene_id, gene_name):
        gid = strip_version(gene_id)
        if gid in id_to_key:
            return id_to_key[gid]
        if gene_name in name_to_key:
            return name_to_key[gene_name]
        return gid  # fall through; classify_position will report unknown_gene

    # ---- classify each RUNFILE row ----
    print("[compute] classifying exon/intron/UTR/intergenic status...", file=sys.stderr)
    feature_types = []
    overlap_transcripts = []
    overlapping_genes_col = []
    position_vs_tss_col = []

    for row in run.itertuples(index=False):
        chrom = str(row.chrom)
        pos = int(row.snp_chromEnd)  # 1-based position, consistent with tss_chromEnd
        gene_key = resolve_gene_key(row.gene_id, row.gene_name)

        ftype, tids = classify_position(gene_key, pos, gene_span, exon_trees, cds_trees, utr_trees)
        feature_types.append(ftype)
        overlap_transcripts.append(tids)

        hits = genes_overlapping(chrom, pos, gene_span)
        overlapping_genes_col.append(";".join(sorted(hits)) if hits else "")

        # strand-aware position relative to TSS
        tss = int(row.tss_chromEnd)
        if pos == tss:
            rel = "at_tss"
        elif row.strand == "+":
            rel = "upstream" if pos < tss else "downstream"
        elif row.strand == "-":
            rel = "upstream" if pos > tss else "downstream"
        else:
            rel = "unknown_strand"
        position_vs_tss_col.append(rel)

    run = run.copy()
    run["feature_type"] = feature_types
    run["overlapping_transcripts"] = overlap_transcripts
    run["genes_overlapping_variant"] = overlapping_genes_col
    run["position_vs_tss"] = position_vs_tss_col

    # ---- merge everything ----
    print("[merge] joining RUNFILE annotations + z-score summary + per-cell-type z-scores...", file=sys.stderr)
    out = run.merge(z_summary, on="SNP", how="left").merge(z, on="SNP", how="left")

    out.to_csv(args.out, index=False)

    print(f"[done] wrote {len(out):,} rows x {len(out.columns)} columns -> {args.out}", file=sys.stderr)
    print("\n[summary] feature_type counts:", file=sys.stderr)
    print(out["feature_type"].value_counts().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()