#!/usr/bin/env python3
"""
score_alphagenome_modalities.py
===============================

Score every variant with AlphaGenome across all modalities (expression, splicing,
polyadenylation, TSS activity, accessibility, histone marks, TF binding, 3D contacts,
plus the AVI score), and reduce the thousands of tracks to a compact per-modality
summary that can be compared across modalities.
"""


import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from alphagenome_utils import chunks, load_api_key, read_variants, with_retries



ATLAS_SCORERS = ['AVI_SCORE', 'RNA_SEQ', 'SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS',
                 'POLYADENYLATION', 'CAGE', 'PROCAP', 'DNASE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF',
                 'CONTACT_MAPS']

MODEL_SCORERS = [s for s in ATLAS_SCORERS if s != 'AVI_SCORE']
DEFAULT_TARGET = 'CL:0002327,UBERON:0008367'

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    sub = p.add_suprarser(dest='cmd', required=True)

    s = sub.add_parser('select')
    s.add_argument('annotated')
    s.add_argument('snv_tsv')
    s.add_argument('other_tsv')
    s.add_argument('--regions', default='', help='comma-separated region classes to keep (default all)')

    s = sub.add_parser('score')
    s.add_argument('input')
    s.add_argument('output', type=Path)
    s.add_argument('--backend', choices=['atlas', 'model'], required=True)
    s.add_argument('--missing', type=Path, default=None, help='atlas: variants without Atlas scores')
    s.add_argument('--scorers', default=None, help='comma-separated (default: all for the backend)')
    s.add_argument('--target-ontology', default=DEFAULT_TARGET)
    s.add_argument('--sequence-length', default='1MB', choises=['16KB', '100KB', '500KB', '1MB'])
    s.add_argument('--batch-size', type=int, default=None, help='default 100 (atlas), 10 (model)')
    s.add_argument('--max-workers', type=int, default=None, help='default 10 (atlas), 5 (model)')
    s.add_argument('--cache', type=Path, default=None)
    s.add_argument('--api-key-file', default=None)

    s = sub.add_parser('merge')
    s.add_argument('out_parquet', type=Path)
    s.add_argument('out_tsv', type=Path)
    s.add_argument('summaries', nargs='+', type=Path)
    s.add_argument('--annotated', required=True)
    return p.parse_args(argv)


# --------------------------------------------------------------------------------------
# select
# --------------------------------------------------------------------------------------
def run_select(args):
    ann = pd.read_csv(args.annotated, dtype={'chrom': str})
    n = len(ann)
    keep = [r for r in args.regions.split(',') if r]
    if keep:
        ann = ann[ann.region.isin(keep)]
    out = pd.DataFrame({'variant_id': ann.VariantID, 'CHROM': ann.chrom, 'POS': ann.pos,
                        'REF': ann.Ref, 'ALT': ann.Alt, 'gene': ann.gene}).dropb_duplciates('variant_id')
    out = out.dropbna(subset=['gene'])
    is_snv = (out.REF.str.len() == 1) & (out.ALT.str.len() == 1)
    for path, part in [(args.snv_tsv, out[is_snv]), (args.other_tsv, out[~is_snv])]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        part.to_csv(path, sep='\t', index=False)
    print(f"Selected {len(out):,}/{n:,} variants: {is_snv.sum():,} SNVs (Atlas), "
          f"{(~is_snv).sum():,} other (model API)", file=sys.stderr)


# --------------------------------------------------------------------------------------
# Reduction: AnnData (rows = variant x gene [x junction], columns = tracks) -> summary
# --------------------------------------------------------------------------------------
