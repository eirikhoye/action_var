import argparse
import csv
import json
import os
import sys
import time
from collections import Counter

import requests

# Script written by biomni ai to convert hgvs formatted custom file to standard vcf file 
# for downstream analysis with tools like AlphaGenome etc.
# TODO: Make it accept a wider file format types 

# --------------------------------------------------------------
# Constants
# --------------------------------------------------------------
VEP_URL = "https://rest.ensembl.org/vep/human/hgvs"
SEQ_URL = "https://rest.ensembl.org/sequence/region/human"
BATCH_SIZE = 200
REQUEST_DELAY = 1

# Convert older RefSeq transcript versions
TX_REMAP = {
    "NM_000059.3": "NM_000059.4",   # BRCA2
    "NM_007294.3": "NM_007294.4",   # BRCA1
    "NM_024675.3": "NM_024675.4",   # PALB2
}

COMP_MAP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
            "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}

VALID_BASES = set("ATGCN")

# GRCh38 chromosome lengths (for VCF contig header lines)
CONTIG_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991089,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


# --------------------------------------------------------------
# Helpers
# --------------------------------------------------------------
def reverse_complement(seq: str) -> str:
    return "".join(COMP_MAP.get(b, "N") for b in reversed(seq))


def parse_allele_string(allele_str: str) -> tuple[str, str]:
    parts = allele_str.split("/")
    ref = parts[0].replace("-", "")
    alt = parts[1].replace("-", "") if len(parts) > 1 else ""
    return ref, alt


def batch_query_vep(hgvs_notations: list[str],
                    batch_size: int = BATCH_SIZE,
                    delay: float = REQUEST_DELAY) -> list[dict]:
    """
    Query Ensembl VEP REST API in batches.
    Returns a list of result dicts (one per successfully resolved HGVS).
    """
    headers = {"Content-Type": "application/json",
               "Accept": "application/json"}
    all_results = []
    n_batches = (len(hgvs_notations) + batch_size - 1) // batch_size

    for i in range(0, len(hgvs_notations), batch_size):
        batch = hgvs_notations[i:i + batch_size]
        batch_num = i // batch_size + 1
        payload = json.dumps({"hgvs_notations": batch})

        for attempt in range(3):
            try:
                resp = requests.post(VEP_URL, headers=headers,
                                     data=payload, timeout=60)
                if resp.status_code == 200:
                    all_results.extend(resp.json())
                    break
                elif resp.status_code == 429:
                    wait = int(resp.headers.get("Retry-After", 5))
                    print(f" Batch {batch_num}/{n_batches}: rate-limited, "
                          f"waiting {wait}s...")
                    time.sleep(wait)
                    continue
                else:
                    print(f" Batch {batch_num}/{n_batches}: HTTP "
                          f"{resp.status_code} - {resp.text[:200]}")
                    break
            except requests.exceptions.RequestException as e:
                print(f" Batch {batch_num}/{n_batches}: error - {e}")
                if attempt < 2:
                    time.sleep(5)
                    continue
                break

        print(f" Batch {batch_num}/{n_batches}: done ({len(batch)} variants)")
        time.sleep(delay)

    return all_results

def fetch_reference_sequences(chrom_ranges: dict[str, list[int]],
                              padding: int = 10) -> dict[str, dict]:
    """
    Fetch reference DNA sequences from Ensembl for the given chromosome ranges.
    Returns {chrom_num: {"start": int, "end": int, "seq": str}}.
    """
    ref_seqs = {}
    for chrom, (s, e) in chrom_ranges.items():
        fetch_start = s - padding
        fetch_end = e + padding
        url = (f"{SEQ_URL}/{chrom}:{fetch_start}-{fetch_end}:1"
               f"?content-type=text/plain")
        for attempt in range(3):
            try:
                resp = requests.get(url, timeout=30)
                if resp.status_code == 200:
                    ref_seqs[chrom] = {
                        "start": fetch_start,
                        "end": fetch_end,
                        "seq": resp.text.strip().upper(),
                    }
                    print(f"  Fetched chr{chrom}: {len(ref_seqs[chrom]['seq'])} bp")
                    break
                elif resp.status_code == 429:
                    time.sleep(int(resp.headers.get("Retry-After", 5)))
                    continue
                else:
                    print(f"  chr{chrom}: HTTP {resp.status_code}")
                    break
            except requests.exceptions.RequestException as e:
                print(f"  chr{chrom}: error - {e}")
                if attempt < 2:
                    time.sleep(3)
        time.sleep(0.5)
    return ref_seqs


def make_ref_base_lookup(ref_seqs: dict[str, dict]):
    """Return a function(pos) → base for a given chrom's reference sequence."""
    def get_ref_base(chrom: str, pos: int) -> str | None:
        r = ref_seqs.get(chrom)
        if r is None:
            return None
        offset = pos - r["start"]
        if 0 <= offset < len(r["seq"]):
            return r["seq"][offset]
        return None
    return get_ref_base


# ── Core conversion ────────────────────────────────────────────────────────

def build_vcf_record(hgvs_id: str, vep_result: dict,
                     gene: str, hgvs_c: str, hgvs_p: str, cls: str,
                     get_ref_base) -> dict | None:
    """
    Convert a single VEP result to a VCF record dict.

    Handles:
      - Strand: VEP returns coding-strand alleles for HGVS input;
        reverse-complement for minus-strand genes to get forward-strand alleles.
      - Indel padding: VCF requires non-empty REF/ALT with a padding base
        for pure insertions and deletions.
    """
    chrom = vep_result["seq_region_name"]
    start = vep_result["start"]
    end = vep_result["end"]
    allele_str = vep_result["allele_string"]
    strand = vep_result.get("strand", 1)

    ref, alt = parse_allele_string(allele_str)

    # Reverse-complement for minus-strand genes (VEP returns coding-strand alleles)
    if strand == -1:
        ref = reverse_complement(ref) if ref else ""
        alt = reverse_complement(alt) if alt else ""

    if ref == "" and alt == "":
        return None

    if ref == "":
        # Insertion: pad with reference base at position `end`
        # (VEP sets start > end for insertions; insertion is between end and end+1)
        pad_pos = end
        pad = get_ref_base(chrom, pad_pos)
        if pad is None:
            return None
        vcf_pos = pad_pos
        vcf_ref = pad
        vcf_alt = pad + alt
    elif alt == "":
        # Deletion: pad with reference base at position start-1
        pad_pos = start - 1
        pad = get_ref_base(chrom, pad_pos)
        if pad is None:
            return None
        vcf_pos = pad_pos
        vcf_ref = pad + ref
        vcf_alt = pad
    else:
        # SNV, MNV, or delins (both alleles non-empty — no padding needed)
        vcf_pos = start
        vcf_ref = ref
        vcf_alt = alt

    # INFO field
    info_parts = [f"GENE={gene}", f"HGVS_C={hgvs_c}"]
    if hgvs_p:
        info_parts.append(f"HGVS_P={hgvs_p}")
    info_parts.append(f"CLASS={cls}")
    info_parts.append(f"TRANSCRIPT={hgvs_id.split(':')[0]}")

    # ID (sanitised HGVS c. notation)
    vid = hgvs_c.replace(">", "_").replace(":", "_")

    return {
        "chrom": f"chr{chrom}",
        "pos": vcf_pos,
        "id": vid,
        "ref": vcf_ref,
        "alt": vcf_alt,
        "info": ";".join(info_parts),
    }


# ── Validation ─────────────────────────────────────────────────────────────

def validate_ref_alleles(records: list[dict], get_ref_base) -> list[dict]:
    """Check that every REF allele matches the reference genome. Returns mismatches."""
    mismatches = []
    for rec in records:
        chrom_num = rec["chrom"].replace("chr", "")
        ref_genome = ""
        for p in range(rec["pos"], rec["pos"] + len(rec["ref"])):
            base = get_ref_base(chrom_num, p)
            if base is None:
                ref_genome = "?"
                break
            ref_genome += base
        if ref_genome != rec["ref"]:
            mismatches.append({
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "ref_vcf": rec["ref"],
                "ref_genome": ref_genome,
                "alt": rec["alt"],
            })
    return mismatches


# ── Writers ────────────────────────────────────────────────────────────────

def write_vcf(records: list[dict], path: str):
    """Write a standard VCFv4.2 file."""
    # Sort by chromosome (numeric) then position
    chrom_nums = {r["chrom"] for r in records}
    chrom_order = {c: int(c.replace("chr", "").replace("X", "23")
                          .replace("Y", "24").replace("M", "25"))
                   for c in chrom_nums}
    records_sorted = sorted(records, key=lambda r: (
        chrom_order.get(r["chrom"], 99), r["pos"], r["ref"], r["alt"]))

    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
        f.write("##source=convert_hgvs_to_vcf.py (via Ensembl VEP REST API)\n")
        f.write("##reference=GRCh38\n")
        for chrom in sorted(chrom_order, key=lambda c: chrom_order[c]):
            length = CONTIG_LENGTHS.get(chrom, "")
            length_str = f",length={length}" if length else ""
            f.write(f"##contig=<ID={chrom}{length_str},assembly=GRCh38>\n")
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
        f.write('##INFO=<ID=HGVS_C,Number=1,Type=String,Description="HGVS coding notation">\n')
        f.write('##INFO=<ID=HGVS_P,Number=1,Type=String,Description="HGVS protein notation">\n')
        f.write('##INFO=<ID=CLASS,Number=1,Type=String,Description="Variant classification">\n')
        f.write('##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="RefSeq transcript accession">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rec in records_sorted:
            info = rec["info"].replace(" ", "_")
            f.write(f"{rec['chrom']}\t{rec['pos']}\t{rec['id']}\t"
                    f"{rec['ref']}\t{rec['alt']}\t.\tPASS\t{info}\n")
    return records_sorted


def write_alphagenome_tsv(records_sorted: list[dict], path: str):
    """Write AlphaGenome batch-scoring TSV (variant_id, CHROM, POS, REF, ALT)."""
    with open(path, "w") as f:
        f.write("variant_id\tCHROM\tPOS\tREF\tALT\n")
        for rec in records_sorted:
            vid = f"{rec['chrom']}_{rec['pos']}_{rec['ref']}_{rec['alt']}_b38"
            f.write(f"{vid}\t{rec['chrom']}\t{rec['pos']}\t"
                    f"{rec['ref']}\t{rec['alt']}\n")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convert HGVS c. notation CSV to VCF + AlphaGenome TSV.")
    parser.add_argument("input_csv", help="Input CSV with Gene, Transcript, "
                        "HGVSc, HGVSp, Class columns")
    parser.add_argument("-o", "--output-prefix", default="variants",
                        help="Output file prefix (default: variants)")
    args = parser.parse_args()

    # 1. Read CSV
    with open(args.input_csv) as f:
        rows = list(csv.DictReader(f))
    print(f"Read {len(rows)} variants from {args.input_csv}")

    # 2. Build HGVS notations, remapping old transcript versions
    hgvs_list = []
    original_hgvs = []  # keep original for mapping back
    for r in rows:
        tx = r["Transcript"]
        c_not = r["HGVSc"]
        orig = f"{tx}:{c_not}"
        new_tx = TX_REMAP.get(tx, tx)
        hgvs_list.append(f"{new_tx}:{c_not}")
        original_hgvs.append(orig)

    # 3. Query Ensembl VEP
    print(f"\nQuerying Ensembl VEP for {len(hgvs_list)} variants...")
    results = batch_query_vep(hgvs_list)

    # Build result map: original HGVS → VEP result
    result_map = {}
    for r in results:
        result_map[r["id"]] = r
    # Map retry results back to original .3 notation
    for orig, new_hgvs in zip(original_hgvs, hgvs_list):
        if orig not in result_map and new_hgvs in result_map:
            result_map[orig] = result_map[new_hgvs]

    missing = [h for h in original_hgvs if h not in result_map]
    if missing:
        # Retry missing with remapped transcripts
        retry_hgvs = []
        retry_orig = []
        for h in missing:
            tx, c_not = h.split(":", 1)
            new_tx = TX_REMAP.get(tx, tx)
            retry_hgvs.append(f"{new_tx}:{c_not}")
            retry_orig.append(h)
        if retry_hgvs:
            print(f"\nRetrying {len(retry_hgvs)} missing variants with remapped transcripts...")
            retry_results = batch_query_vep(retry_hgvs)
            for r, orig in zip(retry_results, retry_orig):
                if "error" not in r:
                    result_map[orig] = r

    still_missing = [h for h in original_hgvs if h not in result_map]
    print(f"\nResolved: {len(result_map)}/{len(original_hgvs)}")
    if still_missing:
        print(f"WARNING: {len(still_missing)} variants could not be resolved:")
        for h in still_missing[:10]:
            print(f"  {h}")

    # 4. Fetch reference sequences for indel padding
    print("\nFetching reference sequences...")
    chrom_ranges = {}
    for h in original_hgvs:
        r = result_map.get(h)
        if r is None:
            continue
        chrom = r["seq_region_name"]
        s, e = r["start"], r["end"]
        if chrom not in chrom_ranges:
            chrom_ranges[chrom] = [s, e]
        else:
            chrom_ranges[chrom][0] = min(chrom_ranges[chrom][0], s)
            chrom_ranges[chrom][1] = max(chrom_ranges[chrom][1], e)

    ref_seqs = fetch_reference_sequences(chrom_ranges)
    get_ref_base = make_ref_base_lookup(ref_seqs)

    # 5. Build VCF records
    records = []
    errors = []
    for i, row in enumerate(rows):
        orig = original_hgvs[i]
        vep = result_map.get(orig)
        if vep is None:
            errors.append(f"{orig}: no VEP result")
            continue
        rec = build_vcf_record(
            orig, vep,
            gene=row["Gene"],
            hgvs_c=row["HGVSc"],
            hgvs_p=row.get("HGVSp", ""),
            cls=row["Class"],
            get_ref_base=get_ref_base,
        )
        if rec is None:
            errors.append(f"{orig}: could not build record")
        else:
            records.append(rec)

    print(f"\nBuilt {len(records)} VCF records ({len(errors)} errors)")

    # 6. Validate REF alleles against reference genome
    mismatches = validate_ref_alleles(records, get_ref_base)
    if mismatches:
        print(f"WARNING: {len(mismatches)} REF alleles do not match reference genome!")
        for m in mismatches[:5]:
            print(f"  {m['chrom']}:{m['pos']} VCF={m['ref_vcf']} genome={m['ref_genome']}")
    else:
        print("All REF alleles verified against GRCh38 reference ✓")

    # 7. Write output files
    vcf_path = f"{args.output_prefix}.vcf"
    tsv_path = f"{args.output_prefix}_alphagenome.tsv"

    records_sorted = write_vcf(records, vcf_path)
    write_alphagenome_tsv(records_sorted, tsv_path)

    print(f"\n── Output files ──")
    print(f"  VCF:              {vcf_path} ({len(records_sorted)} records)")
    print(f"  AlphaGenome TSV:  {tsv_path}")

    # Summary
    gene_counts = Counter()
    for rec in records_sorted:
        gene = rec["info"].split("GENE=")[1].split(";")[0]
        gene_counts[gene] += 1
    print(f"\n── Summary ──")
    for gene, n in sorted(gene_counts.items()):
        print(f"  {gene}: {n} variants")
    print(f"  Total: {len(records_sorted)} variants")


if __name__ == "__main__":
    main()
