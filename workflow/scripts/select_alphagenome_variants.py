#!/usr/bin/env python3
"""
select_alphagenome_variants.py
==============================

Select the variants to score with AlphaGenome from the classify_variant_effects
table, and write them in the input format of the AlphaGenome scripts
(``variant_id, CHROM, POS, REF, ALT, gene``), split by variant type: the
AlphaGenome Atlas only holds SNVs, everything else goes to the model API.

By default all non-coding variants are selected (``coding == False``), which
keeps canonical splice-site, UTR, intronic, promoter and distal variants.

Usage
-----
::

    python select_alphagenome_variants.py EFFECTS_TSV SNV_TSV OTHER_TSV
        [--include noncoding|all] [--exclude-regions "3' UTR,downstream (≤5 kb)"]
"""

import argparse

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('effects_tsv')
    p.add_argument('snv_tsv')
    p.add_argument('other_tsv')
    p.add_argument('--include', choices=['noncoding', 'all'], default='noncoding')
    p.add_argument('--exclude-regions', default='',
                   help='Comma-separated region classes (as in vep_report.py) to leave out')
    args = p.parse_args(argv)

    df = pd.read_csv(args.effects_tsv, sep='\t', dtype={'chrom': str})
    n = len(df)
    if args.include == 'noncoding':
        df = df[~df.coding.astype(bool)]
    excluded = [r for r in args.exclude_regions.split(',') if r]
    if excluded:
        df = df[~df.region.isin(excluded)]
    out = pd.DataFrame({
        'variant_id': df.VariantID, 'CHROM': df.chrom, 'POS': df.pos,
        'REF': df.ref, 'ALT': df.alt, 'gene': df.gene, 'region': df.region,
    }).drop_duplicates('variant_id')
    is_snv = (out.REF.str.len() == 1) & (out.ALT.str.len() == 1)
    out[is_snv].to_csv(args.snv_tsv, sep='\t', index=False)
    out[~is_snv].to_csv(args.other_tsv, sep='\t', index=False)
    print(f"Selected {len(out)}/{n} variants ({args.include}): "
          f"{is_snv.sum()} SNVs, {(~is_snv).sum()} other")


if __name__ == '__main__':
    main()
