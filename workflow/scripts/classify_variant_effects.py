#!/usr/bin/env python3
"""
classify_variant_effects.py
===========================

Apply the variant-effect criteria from Denise's paper to a VEP-annotated VCF,
one row per variant.

Input
-----
``VEP_VCF``
    Output of annotate_vep.py (``--mane``), ideally run with the SIFT/PolyPhen
    options and the plugins in ``vep.plugins`` (see annotate_vep.py). Tools whose
    CSQ fields are absent are reported as missing and their flags left empty
    (NA), so the script also runs on a VCF annotated without plugins.

Method
------
1. One CSQ entry per variant: the MANE Select entry of a target gene (nearest
   if several), else any MANE Select entry, else the first entry. Variants
   without a target-gene MANE entry get region "distal (>5 kb)", as in
   vep_report.py.
2. Per-tool flags, thresholds from the config (defaults = Denise's table):

   ==================  =====================================  ==================
   Question            Tool (CSQ field)                       Default threshold
   ==================  =====================================  ==================
   Loss of function    Consequence                            LoF term present
   Damaging missense   SIFT                                   deleterious
                       PolyPhen (score)                       >= 0.85
                       AlphaMissense (am_class or dbNSFP      likely pathogenic
                       AlphaMissense_pred)
                       REVEL (REVEL_score or REVEL)           >= 0.6
                       MetaRNN (MetaRNN_score)                >= 0.8
   Splicing disrupted  dbscSNV (ada_score, rf_score)          both >= 0.9
                       MaxEntScan (MaxEntScan_diff)           abs >= 2
                       SpliceAI (max SpliceAI_pred_DS_*)      >= 0.5
   Expression altered  ExPectoSC z-scores (optional table)    abs z > 1
   ==================  =====================================  ==================

   Multi-valued fields (``a&b``) use the maximum score / any matching class.
3. Group flags: ``damaging_missense`` needs >= ``missense_min_tools`` tools and
   (by default) a missense consequence; ``splice_disrupted`` needs
   >= ``splice_min_tools`` tools. ``coding`` is True for any coding consequence
   (vep_report.CODING_TERMS) on the chosen entry.

Output
------
``OUT_TSV``   one row per variant: position, gene (target gene; with ``--regions``
              the nearest target gene also for distal variants), transcript, consequence,
              region, coding, raw tool values, per-tool and group flags,
              ``evidence`` (semicolon list of the flags that are True).
``--summary`` optional TSV: per gene, number of variants per group flag.

Usage
-----
::

    python classify_variant_effects.py VEP_VCF OUT_TSV --genes G1,G2,...
        [--thresholds JSON] [--expectosc CSV --expectosc-id-col SNP]
        [--summary TSV]

Author: Eirik Høye
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vep_report import CODING_TERMS, region_class  # noqa: E402  (same CSQ conventions)

DEFAULTS = {
    'lof_consequences': ['transcript_ablation', 'splice_acceptor_variant',
                         'splice_donor_variant', 'stop_gained', 'frameshift_variant',
                         'start_lost'],
    'sift_classes': ['deleterious'],
    'polyphen_min': 0.85,
    'alphamissense_classes': ['likely_pathogenic'],
    'revel_min': 0.6,
    'metarnn_min': 0.8,
    'dbscsnv_min': 0.9,
    'maxentscan_abs_min': 2.0,
    'spliceai_min': 0.5,
    'expectosc_abs_z': 1.0,
    'missense_min_tools': 1,
    'missense_require_consequence': True,
    'splice_min_tools': 1,
}
MISSENSE_TOOLS = ['sift', 'polyphen', 'alphamissense', 'revel', 'metarnn']
SPLICE_TOOLS = ['dbscsnv', 'maxentscan', 'spliceai']
SPLICEAI_FIELDS = ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL',
                   'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']
AM_DBNSFP_CODES = {'P': 'likely_pathogenic', 'B': 'likely_benign', 'A': 'ambiguous'}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('vep_vcf')
    p.add_argument('out_tsv', type=Path)
    p.add_argument('--genes', required=True, help='Comma-separated target genes (HGNC symbols)')
    p.add_argument('--thresholds', default='{}',
                   help='JSON object overriding default thresholds (keys as in DEFAULTS)')
    p.add_argument('--expectosc', default=None,
                   help='Optional ExPectoSC z-score table: one row per variant, one column per cell type')
    p.add_argument('--expectosc-id-col', default='SNP')
    p.add_argument('--expectosc-celltypes', default=None,
                   help='Comma-separated cell-type columns to use (default: all numeric columns)')
    p.add_argument('--regions', default=None,
                   help="vep_report.py --annotated CSV; its gene (nearest target gene for distal "
                        "variants), region and tss_distance replace the CSQ-based ones")
    p.add_argument('--summary', type=Path, default=None)
    return p.parse_args(argv)


# --------------------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------------------
def read_vep_vcf(path):
    """Return (records, csq_fields, csq): records = one row per VCF record,
    csq = one row per CSQ entry, keyed by record index."""
    fields, records, entries = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ'):
                fields = line.split('Format: ')[1].split('"')[0].split('|')
                continue
            if line.startswith('#'):
                continue
            if fields is None:
                raise ValueError(f"{path} has no CSQ header line; is it a VEP VCF?")
            chrom, pos, vid, ref, alt, _, _, info = line.rstrip('\n').split('\t')[:8]
            idx = len(records)
            records.append((vid, chrom, int(pos), ref, alt))
            csq = next((x[4:] for x in info.split(';') if x.startswith('CSQ=')), '')
            for entry in csq.split(',') if csq else []:
                values = entry.split('|')
                if len(values) == len(fields):
                    entries.append([idx] + values)
    records = pd.DataFrame(records, columns=['VariantID', 'chrom', 'pos', 'ref', 'alt'])
    csq = pd.DataFrame(entries, columns=['_rec'] + fields)
    return records, fields, csq


def pick_entries(csq, genes):
    """One CSQ entry per record: target-gene MANE (nearest) > any MANE > first."""
    csq = csq.copy()
    mane = csq.get('MANE_SELECT', pd.Series('', index=csq.index)).fillna('') != ''
    target = csq.SYMBOL.isin(genes) & mane
    csq['_rank'] = np.select([target, mane], [0, 1], 2)
    csq['_dist'] = pd.to_numeric(csq.get('DISTANCE', pd.Series('', index=csq.index))
                                 .replace('', '0'), errors='coerce').fillna(0)
    csq['_order'] = np.arange(len(csq))
    csq = csq.sort_values(['_rec', '_rank', '_dist', '_order'])
    return csq.drop_duplicates('_rec').set_index('_rec')


# --------------------------------------------------------------------------------------
# Value parsing
# --------------------------------------------------------------------------------------
_NUM = re.compile(r'-?\d+(?:\.\d+)?(?:[eE]-?\d+)?')


def max_number(value):
    """Largest number in a CSQ value such as '0.3&0.7', 'probably_damaging(0.99)', '.'."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    nums = [float(x) for part in str(value).split('&')
            for x in _NUM.findall(part.split('(')[-1] if '(' in part else part)]
    return max(nums) if nums else np.nan


def max_abs_number(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    nums = [float(x) for part in str(value).split('&') for x in _NUM.findall(part)]
    return max(nums, key=abs) if nums else np.nan


def classes(value):
    """Set of class labels in e.g. 'deleterious(0.01)&tolerated(0.4)'."""
    if not isinstance(value, str) or value in ('', '.', '-'):
        return set()
    return {part.split('(')[0].strip() for part in value.split('&') if part not in ('', '.', '-')}


def col(df, name):
    return df[name] if name in df.columns else pd.Series(np.nan, index=df.index, dtype=object)


def _na_flags(df):
    return pd.Series(pd.NA, index=df.index, dtype='boolean')


def flag(values, condition):
    """Nullable boolean: NA where the tool has no value."""
    out = pd.Series(pd.NA, index=values.index, dtype='boolean')
    has = values.notna()
    out[has] = condition[has].astype(bool)
    return out


# --------------------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------------------
def classify(records, picked, genes, th):
    df = records.join(picked, how='left')
    out = records.copy()
    ref_len, alt_len = df.ref.str.len(), df.alt.str.len()
    out['variant_type'] = np.select([(ref_len == 1) & (alt_len == 1), ref_len > alt_len,
                                     ref_len < alt_len], ['SNV', 'deletion', 'insertion'],
                                    'MNV / complex')
    out['gene'] = col(df, 'SYMBOL')
    out['transcript'] = col(df, 'Feature')
    out['mane_select'] = col(df, 'MANE_SELECT')
    out['consequence'] = col(df, 'Consequence').fillna('')
    out['impact'] = col(df, 'IMPACT')

    in_target = (df.SYMBOL.isin(genes) & (col(df, 'MANE_SELECT').fillna('') != '')).fillna(False)
    out['region'] = [region_class(c, d or '0')[0] if t else 'distal (>5 kb)'
                     for c, d, t in zip(out.consequence, col(df, 'DISTANCE').fillna('0'), in_target)]
    terms = out.consequence.str.split('&').apply(set)
    out['coding'] = terms.apply(lambda t: bool(t & set(CODING_TERMS)))
    out['missense'] = terms.apply(lambda t: 'missense_variant' in t)

    available = []

    # Loss of function
    lof = set(th['lof_consequences'])
    out['lof'] = terms.apply(lambda t: bool(t & lof)).astype('boolean')

    # Damaging missense
    if 'SIFT' in df:
        available.append('sift')
        out['sift'] = col(df, 'SIFT').replace({'': np.nan})
        out['sift_flag'] = flag(out.sift, out.sift.apply(classes)
                                .apply(lambda c: bool(c & set(th['sift_classes']))))
    if 'PolyPhen' in df:
        available.append('polyphen')
        out['polyphen'] = col(df, 'PolyPhen').apply(max_number)
        out['polyphen_flag'] = flag(out.polyphen, out.polyphen >= th['polyphen_min'])
    if 'am_class' in df or 'AlphaMissense_pred' in df:
        available.append('alphamissense')
        if 'am_class' in df:
            am_cls = col(df, 'am_class').apply(classes)
            out['alphamissense'] = col(df, 'am_pathogenicity').apply(max_number)
        else:
            am_cls = col(df, 'AlphaMissense_pred').apply(classes).apply(
                lambda c: {AM_DBNSFP_CODES.get(x, x) for x in c})
            out['alphamissense'] = col(df, 'AlphaMissense_score').apply(max_number)
        out['alphamissense_class'] = am_cls.apply(lambda c: '&'.join(sorted(c)) or np.nan)
        out['alphamissense_flag'] = flag(out.alphamissense_class, am_cls.apply(
            lambda c: bool(c & set(th['alphamissense_classes']))))
    revel_field = 'REVEL_score' if 'REVEL_score' in df else ('REVEL' if 'REVEL' in df else None)
    if revel_field:
        available.append('revel')
        out['revel'] = col(df, revel_field).apply(max_number)
        out['revel_flag'] = flag(out.revel, out.revel >= th['revel_min'])
    if 'MetaRNN_score' in df:
        available.append('metarnn')
        out['metarnn'] = col(df, 'MetaRNN_score').apply(max_number)
        out['metarnn_flag'] = flag(out.metarnn, out.metarnn >= th['metarnn_min'])

    # Splicing
    if 'ada_score' in df and 'rf_score' in df:
        available.append('dbscsnv')
        out['dbscsnv_ada'] = col(df, 'ada_score').apply(max_number)
        out['dbscsnv_rf'] = col(df, 'rf_score').apply(max_number)
        both = out.dbscsnv_ada.notna() & out.dbscsnv_rf.notna()
        out['dbscsnv_flag'] = flag(out.dbscsnv_ada.where(both),
                                   (out.dbscsnv_ada >= th['dbscsnv_min'])
                                   & (out.dbscsnv_rf >= th['dbscsnv_min']))
    if 'MaxEntScan_diff' in df:
        available.append('maxentscan')
        out['maxentscan_diff'] = col(df, 'MaxEntScan_diff').apply(max_abs_number)
        out['maxentscan_flag'] = flag(out.maxentscan_diff,
                                      out.maxentscan_diff.abs() >= th['maxentscan_abs_min'])
    if any(f in df for f in SPLICEAI_FIELDS):
        available.append('spliceai')
        ds = pd.concat([col(df, f).apply(max_number) for f in SPLICEAI_FIELDS], axis=1)
        out['spliceai_max_ds'] = ds.max(axis=1)
        out['spliceai_flag'] = flag(out.spliceai_max_ds, out.spliceai_max_ds >= th['spliceai_min'])

    # Group flags
    def n_true(tools):
        cols = [f'{t}_flag' for t in tools if f'{t}_flag' in out]
        if not cols:
            return pd.Series(pd.NA, index=out.index, dtype='Int64'), cols
        return out[cols].fillna(False).astype(int).sum(axis=1).astype('Int64'), cols

    out['n_missense_tools'], miss_cols = n_true(MISSENSE_TOOLS)
    dm = out.n_missense_tools >= th['missense_min_tools']
    if th['missense_require_consequence']:
        dm = dm & out.missense
    out['damaging_missense'] = dm.astype('boolean') if miss_cols else _na_flags(out)
    out['n_splice_tools'], splice_cols = n_true(SPLICE_TOOLS)
    out['splice_disrupted'] = ((out.n_splice_tools >= th['splice_min_tools']).astype('boolean')
                               if splice_cols else _na_flags(out))
    return out, available


def add_expectosc(out, path, id_col, celltypes, abs_z):
    z = pd.read_csv(path, sep=None, engine='python')
    if id_col not in z:
        raise ValueError(f"{path} has no column {id_col!r}; set --expectosc-id-col")
    z = z.drop_duplicates(id_col).set_index(id_col)
    cells = celltypes.split(',') if celltypes else list(z.select_dtypes('number').columns)
    z = z[cells]
    abs_max = z.abs().max(axis=1)
    summary = pd.DataFrame({
        'expectosc_max_abs_z': abs_max,
        'expectosc_max_celltype': z.abs().idxmax(axis=1),
        'expectosc_n_celltypes': (z.abs() > abs_z).sum(axis=1),
    })
    out = out.merge(summary, left_on='VariantID', right_index=True, how='left')
    out['expression_altered'] = flag(out.expectosc_max_abs_z, out.expectosc_max_abs_z > abs_z)
    n_hit = out.expectosc_max_abs_z.notna().sum()
    print(f"ExPectoSC: matched {n_hit}/{len(out)} variants on {id_col}", file=sys.stderr)
    return out


GROUP_FLAGS = ['lof', 'damaging_missense', 'splice_disrupted', 'expression_altered']


def main(argv=None):
    args = parse_args(argv)
    th = {**DEFAULTS, **json.loads(args.thresholds)}
    unknown = set(json.loads(args.thresholds)) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"Unknown threshold key(s): {sorted(unknown)}; allowed: {sorted(DEFAULTS)}")
    genes = args.genes.split(',')

    records, fields, csq = read_vep_vcf(args.vep_vcf)
    picked = pick_entries(csq, genes)
    out, available = classify(records, picked, genes, th)
    if args.regions:
        reg = (pd.read_csv(args.regions)[['VariantID', 'gene', 'region', 'tss_distance']]
               .drop_duplicates('VariantID'))
        out = (out.rename(columns={'gene': 'vep_symbol', 'region': 'vep_region'})
                  .merge(reg, on='VariantID', how='left'))
        out['gene'] = out.gene.fillna(out.vep_symbol)
        out['region'] = out.region.fillna(out.vep_region)
    if args.expectosc:
        out = add_expectosc(out, args.expectosc, args.expectosc_id_col,
                            args.expectosc_celltypes, th['expectosc_abs_z'])

    missing = sorted(set(MISSENSE_TOOLS + SPLICE_TOOLS) - set(available))
    print(f"{len(out)} variants; tools found: {available}", file=sys.stderr)
    if missing:
        print(f"Tools without CSQ fields (flags left NA): {missing}. "
              "Enable them under vep: in the config.", file=sys.stderr)

    flags = [f for f in GROUP_FLAGS if f in out]
    hits = pd.DataFrame({f: out[f].astype('boolean').fillna(False).astype(bool) for f in flags})
    out['evidence'] = hits.apply(lambda r: ';'.join(f for f in flags if r[f]), axis=1)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep='\t', index=False)

    if args.summary:
        summary = (out.assign(**hits)
                   .groupby('gene', dropna=False)
                   .agg(variants=('VariantID', 'size'),
                        coding=('coding', 'sum'),
                        **{f: (f, 'sum') for f in flags}))
        summary.to_csv(args.summary, sep='\t')


if __name__ == '__main__':
    main()
