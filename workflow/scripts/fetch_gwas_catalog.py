#!/usr/bin/env python3
"""Fetch GWAS Catalog associations for a trait and write a raw TSV.

Run by the Snakemake rule `download_gwas`, or standalone:

    python fetch_gwas_catalog.py --trait "breast cancer" \
        --p-value-threshold 5e-8 --output out.tsv
"""

# Import libraries
import requests
import pandas as pd
import json
import time


# Using the code examples in official REST API tutorial:
# https://colab.research.google.com/drive/1PDFJAWHN7O6WkZ_HteOkOWeAGnsOMRXR?usp=sharing

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
ENDPOINT = "/v2/studies"


def gwas_api_request(endpoint, params=None):
    """
    Use this function to send GET request for the GWAS Catalog REST API.

    Args:
        endpoint (str): The API endpoint to call (comes after BASE_URL, usually to '/v2/studies/).
        params (dict, optional): A dictionary of query parameters. Default is None.
    
    Returns:
        dict: JSON response from API as python dict, or None if fail.
    """

    full_url = f"{BASE_URL}{endpoint}"

    try:
        response = requests.get(full_url, params=params)

        if response.status_code == 200:
            return response.json()
        else:
            print(f"Error: Received status code {response.status_code}")
            print(f"URL: {response.url}")
            print(f"Response: {response.text}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None

print("Helper function 'gwas_api_request' is defined.")


def get_all_variants_for_trait(trait_name):
    """
    
    """
    variants = set()
    current_page = 0
    total_pages = 1

    print(f"--- Starting search for '{trait_name}")

    # Loop through pages of the API result
    while current_page < total_pages:
        endpoint = "/v2/associations"
        params = {
            "efo_trait": trait_name,
            "size": 200,
            "page": current_page
        }