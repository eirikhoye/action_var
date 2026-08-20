"""
Rules for fetching variant-trait associations from the GWAS Catalog REST API
(https://www.ebi.ac.uk/gwas/rest/api, /v2/associations endpoint).

Include this file from the top-level Snakefile:

    include: "workflow/rules/gwas_catalog.smk"

Entirely driven by the `gwas_catalog:` block in config.yaml -- one output
TSV is produced per trait listed under `gwas_catalog.traits`.
"""
import os
import re

GWAS_CFG = config.get("gwas_catalog", {})
GWAS_TRAITS = GWAS_CFG.get("traits", [])
GWAS_OUTDIR = GWAS_CFG.get("output_dir", "results/gwas_catalog")


def _slugify(text: str) -> str:
    """Turn an EFO trait name into a filesystem/wildcard-safe slug,
    e.g. 'breast carcinoma' -> 'breast_carcinoma'."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


# slug -> original trait name, so the rule can look up the real query string
GWAS_TRAIT_SLUGS = {_slugify(t): t for t in GWAS_TRAITS}

wildcard_constraints:
    slug="|".join(re.escape(s) for s in GWAS_TRAIT_SLUGS) if GWAS_TRAIT_SLUGS else "NOMATCH",


def gwas_catalog_outputs():
    """All expected association TSVs -- use this in `rule all` as a target."""
    return expand(
        os.path.join(GWAS_OUTDIR, "{slug}.associations.tsv"),
        slug=GWAS_TRAIT_SLUGS.keys(),
    )


rule fetch_gwas_associations:
    """Fetch every association for one EFO trait from the GWAS Catalog REST API."""
    output:
        tsv=os.path.join(GWAS_OUTDIR, "{slug}.associations.tsv"),
    params:
        trait=lambda wc: GWAS_TRAIT_SLUGS[wc.slug],
        base_url=GWAS_CFG.get("base_url", "https://www.ebi.ac.uk/gwas/rest/api"),
        page_size=GWAS_CFG.get("page_size", 200),
        p_value_threshold=GWAS_CFG.get("p_value_threshold", None),
        sleep=GWAS_CFG.get("request_sleep_seconds", 0.1),
    log:
        "logs/gwas_catalog/{slug}.log",
    retries: 2
    script:
        "../scripts/fetch_gwas_associations.py"