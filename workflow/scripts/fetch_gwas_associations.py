"""
Fetch all GWAS Catalog associations for a given EFO trait via the
GWAS Catalog REST API (/v2/associations) and write them to a TSV file.

Called by workflow/rules/gwas_catalog.smk via Snakemake's `script:` directive.
Expects a `snakemake` object (injected automatically) with:
    params.trait, params.base_url, params.page_size,
    params.p_value_threshold, params.sleep
    output.tsv
    log[0]

Modified from the code examples in official REST API tutorial:
https://colab.research.google.com/drive/1PDFJAWHN7O6WkZ_HteOkOWeAGnsOMRXR?usp=sharing

"""
from __future__ import annotations

import csv
import logging
import time 
import argparse
import re
import requests


#-------------------------------------------------------------------------------
# ARGS
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
                    description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('trait', type=str, help='The trait that is to find associated variants')
parser.add_argument('--base_url', type=str, default='https://www.ebi.ac.uk/gwas/rest/api',
                     help='API base URL (default: %(default)s)')
parser.add_argument('--page_size', type=int, default=200,
                    help="Size per page, larger for more efficiency (default %(default)s)")
parser.add_argument('--p_value_threshold', type=float, default=5.0e-8,
                    help='Minimal p-value for association with trait (default: %(default)s)')
parser.add_argument('--sleep_s', type=float, default=0.1,
                    help='Wait intervals between requests, (default %(default)s')
parser.add_argument('--out_path', type=str, default=None,
                    help='Path to output file (default: <trait>_variants.tsv)')
parser.add_argument('--log_path', type=str, default=None,
                    help='Path to logfile (default: <trait>_variants.log)')

args = parser.parse_args()

trait = args.trait
# 'breast cancer' -> 'breast_cancer', safe for use in a filename
trait_slug = re.sub(r"[^a-z0-9]+", "_", trait.lower()).strip("_")

out_path = args.out_path or f"{trait_slug}_variants.tsv"
log_path = args.log_path or f"{trait_slug}_variants.log"

base_url, page_size, p_value_threshold, sleep_s = (
    args.base_url, args.page_size, args.p_value_threshold, 
    args.sleep_s
)

logging.basicConfig(
    filename=log_path,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)
    

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------
def gwas_api_request(endpoint: str, params: dict) -> dict | None:
    url = f"{base_url}{endpoint}"
    try:
        resp = requests.get(url, params=params, timeout=60)
    except requests.exceptions.RequestException as exc:
        logger.error("Network error calling %s: %s", url, exc)
        return None

    if resp.status_code != 200:
        logger.error("HTTP %s from %s | body: %s", resp.status_code, resp.url, resp.text[:500])
        return None

    return resp.json()


def fetch_all_associations(trait_name:str) -> list[list]:
    """Paginate through /v2/associations for a given EFO trait, sorted by p-value."""
    all_associations: list[dict] = []
    current_page = 0
    total_pages = 1

    logger.info("Starting fetch for trait '%s'", trait_name)

    while current_page < total_pages:
        data = gwas_api_request(
            "/v2/associations",
            {
                "efo_trait": trait_name,
                "sort": "p_value",
                "direction": "asc",
                "size": page_size,
                "page": current_page,
            },
        )

        if not data or "_embedded" not in data:
            logger.warning("No data / error on page %d. Stopping.", current_page)
            break

        associations = data["_embedded"]["associations"]
        all_associations.extend(associations)

        if current_page == 0:
            total_pages = data["page"]["totalPages"]
            logger.info(
                "Found %d total associations across %d pages",
                data["page"]["totalElements"],
                total_pages,
            )

        current_page += 1
        time.sleep(sleep_s)

    logger.info("Finished fetch for '%s': %d associations", trait_name, len(all_associations))
    return all_associations

def parse_row(assoc: dict) -> dict:
    """Flatten one association object into a TSV row, matching the fields
    used in the GWAS Catalog REST API sample notebook."""
    rs_id, risk_allele = None, None
    snp_effect_allele = assoc.get("snp_effect_allele") or []
    if snp_effect_allele:
        risk_allele_str = snp_effect_allele[0]
        if "-" in risk_allele_str:
            rs_id, risk_allele = risk_allele_str.split("-", 1)
        else:
            rs_id = risk_allele_str

    mapped_genes = assoc.get("mapped_genes") or []
    mapped_genes_str = ", ".join(mapped_genes) if isinstance(mapped_genes, list) else str(mapped_genes)

    efo_traits = assoc.get("efo_traits") or []
    efo_trait_names = ", ".join(t.get("efo_trait", "") for t in efo_traits)

    return {
        "variant_rsID": rs_id,
        "risk_allele_base": risk_allele,
        "p_value": assoc.get("p_value"),
        "or_value": assoc.get("or_value"),
        "risk_frequency": assoc.get("risk_frequency"),
        "mapped_genes": mapped_genes_str,
        "efo_traits": efo_trait_names,
        "accession_id": assoc.get("accession_id"),
        "pubmed_id": assoc.get("pubmed_id"),
        "first_author": assoc.get("first_author"),
    }


def main():
    associations = fetch_all_associations(trait)
    rows = [parse_row(a) for a in associations]

    if p_value_threshold is not None:
        before = len(rows)
        rows = [
            r for r in rows
            if r["p_value"] is not None and r["p_value"] <= float(p_value_threshold)
        ]
        logger.info(
            "Applied p-value filter (<= %s): %d -> %d rows",
            p_value_threshold, before, len(rows),
        )

    fieldnames = [
        "variant_rsID", "risk_allele_base", "p_value", "or_value",
        "risk_frequency", "mapped_genes", "efo_traits",
        "accession_id", "pubmed_id", "first_author",
    ]

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote %d rows to %s", len(rows), out_path)


if __name__ == "__main__":
    main()