#!/usr/bin/env python3
"""
Fetch GWAS Catalog associations for a trait and write a raw TSV.

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


def get_all_variants_for_trait(trait_name):
    """
    
    """
    variants = set()
    current_page = 0
    total_pages = 1

    print(f"--- Starting search for '{trait_name}' ---")

    while current_page < total_pages:
        endpoint = "/v2/associations"
        params = {
            "efo_trait": trait_name,
            "size": 200,
            "page": current_page
        }

        data = gwas_api_request(endpoint, params)

        if data and '_embedded' in data:
            associations_list = data['_embedded']['associations']

            if current_page == 0:
                total_pages = data['page']['totalPages']

            for association in associations_list:
                if 'snp_effect_allele' in association and association['snp_effect_allele']:
                    risk_allele_str = association['snp_effect_allele'][0]
                    rs_id = risk_allele_str.split('-')[0]
                    if rs_id.startswith('rs'):
                        variants.add(rs_id)

            current_page += 1
            time.sleep(0.1)
        else:
            print("No more data or an error occurred. Stopping.")
            break

    print(f"--- Finished fetching for '{trait_name}'. Found {len(variants)} total unique variants. ---\n")       
    return variants


def get_all_associations_for_trait(trait_name):
    """
    
    """
    all_associations = []
    current_page = 0
    total_pages = 1

    print(f"--- Starting search for all associations related to '{trait_name}' ---")

    while current_page < total_pages:
        endpoint = "/v2/associations"
        params = {
            "efo_trait": trait_name,
            "sort" :"p_value",
            "direction": "asc",
            "size": 40,
            "page": current_page
        }

        data = gwas_api_request(endpoint, params)

        if data and '_embedded' in data:
            associations_list = data['_embedded']['associations']
            all_associations.extend(associations_list)

            if current_page == 0:
                total_pages = data['page']['totalPages']

            current_page += 1
            time.sleep(0.1)
        else:
            print("No more data or an error occurred. Stopping.")
            break

    print(f"--- Finished fetching. Found {len(all_associations)} total associations. ---\n")
    return all_associations


def get_all_genes_for_trait(trait_name):
    all_genes = set()
    current_page = 0
    total_pages = 1

    while current_page < total_pages:
        endpoint = "/v2/associations"
        params = {
            "efo_trait": trait_name,
            "size": 200,
            "page": current_page
        }

        associations_data = gwas_api_request(endpoint, params)

        if associations_data and '_embedded' in associations_data:
            associations_list = associations_data['_embedded']['associations']

            if current_page == 0:
                total_pages = associations_data['page']['totalPages']

            for association in associations_list:
                if 'mapped_genes' in association and association['mapped_genes']:
                    for gene_string in association['mapped_genes']:
                        genes = gene_string.split(',')
                        for gene in genes:
                            if gene.strip():
                                all_genes.add(gene.strip())

            current_page += 1
        else:
            print("No more data found or an error occurred. Stopping.")
            break

        return all_genes


def get_all_associations_for_accession_id(accession_id):
    """
    
    """
    all_associations = []
    current_page = 0
    total_pages = 1

    print(f"--- Starting search for all associations related to '{accession_id}' ---")

    while current_page < total_pages:
        endpoint = "/v2/associations"
        params = {
            "accession_id": accession_id,
            "sort": "p_value",
            "direction": "asc",
            "size": 200,
            "page": current_page
        }

        data = gwas_api_request(endpoint, params)

        if data and '_embedded' in data:
            associations_list = data['_embedded']['associations']
            all_associations.extend(associations_list)

            if current_page == 0:
                total_pages = data['page']['totalPages']
                print(f"Found {data['page']['totalElements']} total associations across {total_pages} pages.")

            print(f"Page {current_page + 1}/{total_pages} processed. Collected {len(all_associations)} associations so far.")
            current_page += 1
            time.sleep(0.1)
        else:
            print("No more data or an error occurred. Stopping")
            break

    print(f"--- Finished fetching. Found {len(all_associations)} total associations. ---\n")
    return all_associations


