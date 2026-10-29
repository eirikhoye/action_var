#!/usr/bin/env python3
"""
merge_sei_scores.py
===================

Join Sei sequence-class variant effects back onto the variant table.

Input
-----
``SCORES_NPY``   ``<prefix>.sequence_class_scores.npy`` from
                 ``2_varianteffect_sc_score.py --no-tsv`` (alt - ref, histone-
                 normalised), rows in the order of ...
``ROW_LABELS``   ``<prefix>_row_labels.txt`` (Selene: chrom, pos, name, ref, alt,
                 strand, ref_match, contains_unk)
``SEQCLASS``     ``model/seqclass.names`` from the Sei framework (one class per line,
                 e.g. ``TF3 FOXA1 / AR / ESR1``)
``EFFECTS_TSV``  classify_variant_effects output (``VariantID`` column)

Output
------
``OUT_TSV``: the effects table (all rows) plus, per variant:

* ``sei_ref_match``, ``sei_contains_unk``: Selene QC; a False ref_match means the
  VCF REF did not match the FASTA, so treat those scores with care
* ``sei_max_abs``, ``sei_top_class``, ``sei_top_score``: strongest class effect
* ``sei_<group>_max``: signed score of the largest-magnitude class per group
  (P, E, CTCF, TF, PC, HET, TN, L)
* ``sei_<label>`` for every sequence class (e.g. ``sei_TF3``, ``sei_E6``)

``--long`` optionally writes a long table (VariantID, class, name, group, score).

Usage
-----
::

    python merge_sei_scores.py SCORES_NPY ROW_LABELS SEQCLASS EFFECTS_TSV OUT_TSV [--long TSV]
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd


def read_classes(path, n_columns):
    names = [line.strip() for line in open(path) if line.strip()]
    labels = [n.split(' ', 1)[0] for n in names]
    if n_columns > len(names):   # unnamed extra classes, if any
        extra = range(len(names), n_columns)
        names += [f'SC{i} unnamed' for i in extra]
        labels += [f'SC{i}' for i in extra]
    return pd.DataFrame({'label': labels[:n_columns], 'name': names[:n_columns],
                         'group': [re.match(r'[A-Z]+', lab).group(0) for lab in labels[:n_columns]]})


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('scores_npy')
    p.add_argument('row_labels')
    p.add_argument('seqclass_names')
    p.add_argument('effects_tsv')
    p.add_argument('out_tsv')
    p.add_argument('--long', default=None)
    args = p.parse_args(argv)

    scores = np.load(args.scores_npy)
    rows = pd.read_csv(args.row_labels, sep='\t', dtype={'chrom': str})
    if len(rows) != scores.shape[0]:
        raise ValueError(f"{args.row_labels} has {len(rows)} rows but {args.scores_npy} "
                         f"has {scores.shape[0]}")
    classes = read_classes(args.seqclass_names, scores.shape[1])

    sc = pd.DataFrame(scores, columns='sei_' + classes.label)
    sc.insert(0, 'VariantID', rows['name'].astype(str).values)
    sc['sei_ref_match'] = rows['ref_match'].values
    sc['sei_contains_unk'] = rows['contains_unk'].values
    absval = np.abs(scores)
    top = absval.argmax(axis=1)
    sc['sei_max_abs'] = absval.max(axis=1)
    sc['sei_top_class'] = classes.name.values[top]
    sc['sei_top_score'] = scores[np.arange(len(scores)), top]
    for group, idx in classes.groupby('group').groups.items():
        sub = scores[:, idx]
        sc[f'sei_{group}_max'] = sub[np.arange(len(sub)), np.abs(sub).argmax(axis=1)]
    sc = sc.drop_duplicates('VariantID')

    effects = pd.read_csv(args.effects_tsv, sep='\t')
    summary_cols = ['sei_ref_match', 'sei_contains_unk', 'sei_max_abs', 'sei_top_class',
                    'sei_top_score'] + sorted(c for c in sc if c.endswith('_max') and c != 'sei_max_abs')
    class_cols = list('sei_' + classes.label)
    out = effects.merge(sc[['VariantID'] + summary_cols + class_cols], on='VariantID', how='left')
    out.to_csv(args.out_tsv, sep='\t', index=False)
    n_scored = out.sei_max_abs.notna().sum()
    n_bad = (sc.sei_ref_match.astype(str) == 'False').sum()
    print(f"Sei scores for {n_scored}/{len(out)} variants; {n_bad} with REF mismatch",
          file=sys.stderr)

    if args.long:
        long = sc.melt(id_vars='VariantID', value_vars=class_cols, var_name='label',
                       value_name='score')
        long['label'] = long.label.str.removeprefix('sei_')
        long = long.merge(classes, on='label')
        long.to_csv(args.long, sep='\t', index=False)


if __name__ == '__main__':
    main()
