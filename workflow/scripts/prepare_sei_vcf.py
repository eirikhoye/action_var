#!/usr/bin/env python3
"""
prepare_sei_vcf.py
==================

Write the minimal VCF that Sei (via Selene) expects, from the
classify_variant_effects table: ``chr``-prefixed CHROM, POS, ID = VariantID,
REF, ALT, then ``. . .``. Sei's bundled hg38 FASTA uses UCSC names, so the
chromosome naming must match.

* ``--include noncoding`` (default) keeps ``coding == False`` variants, as for
  AlphaGenome; ``all`` keeps every variant.
* ``--snv-only`` drops indels. Otherwise indels are kept; Selene centres the
  4,096 bp window on the variant and reports REF mismatches in ``ref_match``.
* Alleles with characters other than ACGT, and duplicate IDs, are dropped.

Usage
-----
::

    python prepare_sei_vcf.py EFFECTS_TSV OUT_VCF [--include noncoding|all] [--snv-only]
"""

import argparse
import sys

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('effects_tsv')
    p.add_argument('out_vcf')
    p.add_argument('--include', choices=['noncoding', 'all'], default='noncoding')
    p.add_argument('--snv-only', action='store_true')
    args = p.parse_args(argv)

    df = pd.read_csv(args.effects_tsv, sep='\t', dtype={'chrom': str})
    n = len(df)
    if args.include == 'noncoding':
        df = df[~df.coding.astype(bool)]
    acgt = df.ref.str.fullmatch('[ACGT]+') & df.alt.str.fullmatch('[ACGT]+')
    if (~acgt).any():
        print(f"Dropping {(~acgt).sum()} variants with non-ACGT alleles", file=sys.stderr)
    df = df[acgt]
    if args.snv_only:
        df = df[(df.ref.str.len() == 1) & (df.alt.str.len() == 1)]
    df = df.drop_duplicates('VariantID').copy()
    df['chrom'] = 'chr' + df.chrom.str.removeprefix('chr')
    df = df.sort_values(['chrom', 'pos'])

    with open(args.out_vcf, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for r in df.itertuples():
            fh.write(f"{r.chrom}\t{r.pos}\t{r.VariantID}\t{r.ref}\t{r.alt}\t.\t.\t.\n")
    print(f"Wrote {len(df)}/{n} variants ({args.include}) to {args.out_vcf}", file=sys.stderr)


if __name__ == '__main__':
    main()
