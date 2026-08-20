#!/usr/bin/env python3
import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ENSEMBL_SERVER = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}
BATCH_SIZE = 200          # Ensembl POST batch limit is generous, but keep requests small
REQUEST_PAUSE = 0.34      # ~3 req/s to stay well under the anonymous rate limit


# ---------------------------------------------------------------------------
# Step 1-2: load, dedup, rescue, split
# ---------------------------------------------------------------------------
def clean_and_dedup(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df["or_value"] = pd.to_numeric(df["or_value"], errors="coerce")
    df["is_coord_id"] = ~df["variant_rsID"].str.startswith("rs")

    agg = (
        df.sort_values("p_value")
        .groupby("variant_rsID", as_index=False)
        .agg(
            risk_allele_base=("risk_allele_base", lambda s: next((x for x in s if x != "?"), "?")),
            best_p_value=("p_value", "min"),
            or_value=("or_value", "first"),
            risk_frequency=("risk_frequency", lambda s: next((x for x in s if x != "NR"), "NR")),
            mapped_genes=("mapped_genes", "first"),
            efo_traits=("efo_traits", lambda s: "; ".join(sorted(set(s)))),
            accession_ids=("accession_id", lambda s: ";".join(sorted(set(s)))),
            n_studies=("accession_id", "nunique"),
            is_coord_id=("is_coord_id", "first"),
        )
    )
    return agg


# ---------------------------------------------------------------------------
# Step 3: resolve rsIDs -> SPDI via Ensembl variant_recoder
# ---------------------------------------------------------------------------
def resolve_rsids(rsids, assembly="GRCh38"):
    """Returns dict: rsid -> list of (chrom, pos_1based, ref, alt) candidate alleles."""
    server = ENSEMBL_SERVER[assembly]
    url = f"{server}/variant_recoder/human"
    resolved = {}
    rsids = list(rsids)

    for i in range(0, len(rsids), BATCH_SIZE):
        batch = rsids[i:i + BATCH_SIZE]
        try:
            r = requests.post(
                url,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                json={"ids": batch},
                timeout=60,
            )
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"  [warn] batch {i}-{i+len(batch)} failed: {e}", file=sys.stderr)
            time.sleep(2)
            continue

        for entry in r.json():
            # each entry is {input_rsid: {...spdi results...}} OR has 'input' key depending on version
            for rsid, payload in entry.items():
                if rsid == "warnings":
                    continue
                spdi_list = payload.get("spdi", [])
                alleles = []
                for spdi in spdi_list:
                    # SPDI format: "NC_000010.11:123279685:A:G"  (0-based interbase position)
                    try:
                        seq_id, pos0, ref, alt = spdi.split(":")
                        chrom = _refseq_to_chrom(seq_id)
                        pos1 = int(pos0) + 1  # convert 0-based interbase -> 1-based VCF-style
                        if ref:  # skip pure insertions with empty ref for the simple SNP case
                            alleles.append((chrom, pos1, ref, alt))
                    except ValueError:
                        continue
                if alleles:
                    resolved[rsid] = alleles
        time.sleep(REQUEST_PAUSE)

    return resolved


def _refseq_to_chrom(refseq_id):
    """Map RefSeq chromosome accessions (NC_000001.x ... NC_000024.x) to chr1..chrY."""
    mapping = {f"NC_{i:06d}": str(i) for i in range(1, 23)}
    mapping["NC_000023"] = "X"
    mapping["NC_000024"] = "Y"
    prefix = refseq_id.split(".")[0]
    return "chr" + mapping.get(prefix, prefix)


# ---------------------------------------------------------------------------
# Step 5: reference base lookup for coordinate-only entries
# ---------------------------------------------------------------------------
def fetch_ref_base(chrom, pos, assembly="GRCh38"):
    server = ENSEMBL_SERVER[assembly]
    region = f"{chrom.replace('chr', '')}:{pos}-{pos}"
    url = f"{server}/sequence/region/human/{region}"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
        r.raise_for_status()
        return r.json().get("seq")
    except requests.RequestException as e:
        print(f"  [warn] ref lookup failed for {chrom}:{pos}: {e}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Step 4/6: match risk allele, build final table
# ---------------------------------------------------------------------------
def build_resolved_table(agg: pd.DataFrame, assembly="GRCh38"):
    rs_rows = agg[~agg["is_coord_id"]].copy()
    coord_rows = agg[agg["is_coord_id"]].copy()

    print(f"Resolving {rs_rows['variant_rsID'].nunique()} rsIDs via Ensembl variant_recoder ...")
    spdi_map = resolve_rsids(rs_rows["variant_rsID"], assembly=assembly)

    resolved_records = []
    unresolved_records = []

    for _, row in rs_rows.iterrows():
        rsid = row["variant_rsID"]
        risk_allele = row["risk_allele_base"]
        candidates = spdi_map.get(rsid)

        if not candidates:
            unresolved_records.append({**row.to_dict(), "reason": "no Ensembl coordinate mapping"})
            continue

        if risk_allele == "?":
            unresolved_records.append({**row.to_dict(), "reason": "risk allele unknown ('?') in GWAS Catalog export",
                                        "candidate_alleles": candidates})
            continue

        match = next((c for c in candidates if c[3] == risk_allele), None)
        if match is None:
            unresolved_records.append({**row.to_dict(),
                                        "reason": f"risk allele '{risk_allele}' not among Ensembl alts {[c[3] for c in candidates]}",
                                        "candidate_alleles": candidates})
            continue

        chrom, pos, ref, alt = match
        resolved_records.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "variant_id": rsid, **_shared_fields(row),
        })

    for _, row in coord_rows.iterrows():
        chrom, pos = row["variant_rsID"].split(":")
        pos = int(pos)
        risk_allele = row["risk_allele_base"]
        if risk_allele == "?":
            unresolved_records.append({**row.to_dict(), "reason": "risk allele unknown ('?')"})
            continue
        ref = fetch_ref_base(chrom, pos, assembly=assembly)
        time.sleep(REQUEST_PAUSE)
        if ref is None:
            unresolved_records.append({**row.to_dict(), "reason": "reference base lookup failed"})
            continue
        if ref == risk_allele:
            unresolved_records.append({**row.to_dict(),
                                        "reason": "risk allele == reference base (can't build ALT without a second known allele)"})
            continue
        resolved_records.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": risk_allele,
            "variant_id": row["variant_rsID"], **_shared_fields(row),
        })

    resolved_df = pd.DataFrame(resolved_records)
    unresolved_df = pd.DataFrame(unresolved_records)
    return resolved_df, unresolved_df


def _shared_fields(row):
    return {
        "gene": row["mapped_genes"],
        "trait": row["efo_traits"],
        "p_value": row["best_p_value"],
        "or_value": row["or_value"],
        "risk_frequency": row["risk_frequency"],
        "n_studies": row["n_studies"],
        "accession_ids": row["accession_ids"],
    }


def write_vcf(resolved_df: pd.DataFrame, path: Path):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=resolve_gwas_variants.py\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, r in resolved_df.sort_values(["chrom", "pos"]).iterrows():
            info = f"GENE={r['gene']};TRAIT={r['trait']};OR={r['or_value']};P={r['p_value']}"
            fh.write(f"{r['chrom']}\t{r['pos']}\t{r['variant_id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t{info}\n")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_tsv")
    ap.add_argument("--outdir", default="resolved")
    ap.add_argument("--assembly", choices=["GRCh38", "GRCh37"], default="GRCh38")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    agg = clean_and_dedup(args.input_tsv)
    print(f"Loaded and deduped to {len(agg)} unique variants.")

    resolved_df, unresolved_df = build_resolved_table(agg, assembly=args.assembly)

    resolved_path = outdir / "variants_resolved.tsv"
    unresolved_path = outdir / "variants_unresolved.tsv"
    vcf_path = outdir / "variants_resolved.vcf"

    resolved_df.to_csv(resolved_path, sep="\t", index=False)
    unresolved_df.to_csv(unresolved_path, sep="\t", index=False)
    if len(resolved_df):
        write_vcf(resolved_df, vcf_path)

    print(f"\nResolved:   {len(resolved_df):4d}  -> {resolved_path}")
    print(f"Unresolved: {len(unresolved_df):4d}  -> {unresolved_path}")
    if len(resolved_df):
        print(f"VCF:                -> {vcf_path}")


if __name__ == "__main__":
    main()










