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

import requests


#-------------------------------------------------------------------------------
# ARGS
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Fetch all GWAS Catalog associations for a given EFO trait")
parser.add_argument('trait', type=str, help='The trait that is to find associated variants')
parser.add_argument('base_url', type=str, help='URL, typically: https://www.ebi.ac.uk/gwas/rest/api')
parser.add_argument('page_size', type=int, help="Size per page, larger for more efficiency")
parser.add_argument('p_value_threshold', type=float, help='Minimal p-value for association with trait')
parser.add_argument('sleep_s', type=float, help='Wait intervals between requests, typical 0.1')
parser.add_argument('out_path', type=str, help='Path to output file')
parser.add_argument('log_path', type=str, help='Path to logfile')
args = parser.parse_args()
trait, base_url, page_size, p_value_threshold, sleep_s, out_path, log_path = args.trait, args.base_url, args.page_size, args.p_value_threshold, args.sleep_s, args.out_path, args.log_path

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



#---------------- PLACEHOLDER REOMVE ----------------
endpoint = "/v2/associations"
params = {
    "efo_trait": trait,
    "size": 10
}
#---------------- PLACEHOLDER REOMVE ----------------



def main():
    print(gwas_api_request(endpoint, params))

if __name__ == "__main__":
    main()