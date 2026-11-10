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
    sub = p.add_subparsers(dest='cmd', required=True)

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
    s.add_argument('--sequence-length', default='1MB', choices=['16KB', '100KB', '500KB', '1MB'])
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
                        'REF': ann.Ref, 'ALT': ann.Alt, 'gene': ann.gene}).drop_duplicates('variant_id')
    out = out.dropna(subset=['gene'])
    is_snv = (out.REF.str.len() == 1) & (out.ALT.str.len() == 1)
    for path, part in [(args.snv_tsv, out[is_snv]), (args.other_tsv, out[~is_snv])]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        part.to_csv(path, sep='\t', index=False)
    print(f"Selected {len(out):,}/{n:,} variants: {is_snv.sum():,} SNVs (Atlas), "
          f"{(~is_snv).sum():,} other (model API)", file=sys.stderr)


# --------------------------------------------------------------------------------------
# Reduction: AnnData (rows = variant x gene [x junction], columns = tracks) -> summary
# --------------------------------------------------------------------------------------
def track_labels(var):
    name = var['name'].astype(str)
    if 'biosample_name' in var:
        bio = var['biosample_name'].astype(str)
        name = np.where(var['biosample_name'].notna(), bio + ' | ' + name, name)
    return np.asarray(name, dtype=object)


def reduce_scorer(scorer, ad, row_keys, row_genes, targets):
    """Per variant: strongest target-biosample track and strongest track overall.

    row_keys: variant key of each AnnData row; row_genes: the variant's target gene for
    each row (rows for other genes are dropped). Returns a DataFrame indexed by key.
    """
    X = np.asarray(ad.X, dtype=float)
    Q = (np.asarray(ad.layers['quantiles'], dtype=float) if 'quantiles' in ad.layers
         else np.full_like(X, np.nan))
    keep = np.ones(ad.n_obs, dtype=bool)
    if 'gene_name' in ad.obs:
        keep = (ad.obs['gene_name'].astype(str).values == np.asarray(row_genes, dtype=str))
    if not keep.any():
        return pd.DataFrame()
    X, Q, keys = X[keep], Q[keep], np.asarray(row_keys)[keep]
    obs = ad.obs[keep]
    var = ad.var
    labels = track_labels(var)
    has_bio = 'ontology_curie' in var and var['ontology_curie'].notna().any()
    target_cols = var['ontology_curie'].isin(targets).values if has_bio else np.ones(X.shape[1], bool)

    out = {}
    for suffix, cols in [('', target_cols), ('_all', np.ones(X.shape[1], bool))]:
        if suffix == '_all' and not has_bio:
            continue   # identical to the target value
        if not cols.any():
            continue
        x, q, lab = X[:, cols], Q[:, cols], labels[cols]
        # rank by |quantile|, ties (and missing quantiles) by |raw|
        rank = np.nan_to_num(np.abs(q), nan=-1.0) + 1e-6 * np.tanh(np.abs(np.nan_to_num(x)))
        j = rank.argmax(axis=1)
        i = np.arange(len(j))
        df = pd.DataFrame({'key': keys, 'rank': rank[i, j], 'q': q[i, j], 'raw': x[i, j],
                           'track': lab[j]})
        if 'junction_Start' in obs:
            df['junction'] = (obs['junction_Start'].astype('Int64').astype(str) + '-'
                              + obs['junction_End'].astype('Int64').astype(str)).values
        best = df.sort_values('rank', ascending=False).drop_duplicates('key').set_index('key')
        out[f'{scorer}|q{suffix}'] = best.q
        out[f'{scorer}|raw{suffix}'] = best.raw
        if suffix == '_all':
            out[f'{scorer}|track_all'] = best.track
        if 'junction' in best and suffix == '':
            out[f'{scorer}|junction'] = best.junction
    return pd.DataFrame(out)


def summarise(per_scorer, batch, backend):
    """Join per-scorer summaries (indexed by variant key) onto the batch's variants."""
    from alphagenome.data import genome
    base = batch[['variant_id', 'gene', 'CHROM', 'POS', 'REF', 'ALT']].copy()
    base['key'] = [str(genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT)) for r in batch.itertuples()]
    frames = [f for f in per_scorer if not f.empty]
    if frames:
        base = base.join(pd.concat(frames, axis=1), on='key')
    base['backend'] = backend
    return base.drop(columns='key')    


# --------------------------------------------------------------------------------------
# score: Atlas
# --------------------------------------------------------------------------------------
def atlas_batch(client, batch, gene, scorers, targets, workers):
    from alphagenome.data import genome
    variants = [genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT, name=r.variant_id)
                for r in batch.itertuples()]

    def query(vs):
        return client.query_variants(vs, requested_scorers=scorers, gene_names=[gene],
                                     progress_bar=False, max_wlrkers=workers)

    try:
        results = [with_retries(lambda: query(variants), f'{gene} batch', attempts=3)]
    except ValueError: # Check if one variant is not in Atlas
        results = []
        for v in variants:
            try:
                results.append(query([v]))
            except ValueError as err:
                print(f"not in Atlas: {v} ({str(err)[:100]})", file=sys.stderr)
    per_scorer = []
    for result in results:
        for scorer, ad in results.items():
            if ad.n_obs == 0 or 'variant' not in ad.obs:
                continue
            keys = ad.obs['variant'].astype(str).values
            per_scorer.append(reduce_scorer(scorer, ad, keys, [gene] * ad.n_obs, targets))
    # several results (one-by-one fallback): stack rows per scorer before joining
    merged = {}
    for f in per_scorer:
        for col in f:
            merged.setdefault(col, []).append(f[col])
    frame = pd.DataFrame({c: pd.concat(v)[lambda s: ~s.index.duplicated()] for c, v in merged.items()})
    return summarise([frame], batch, 'atlas')


def atlas_reason(client, r):
    """Why the Atlas has no scores for a variant. A reference mismatch means the variant's
    REF is not the hg38 base (e.g. a liftover problem): the model API would score it anyway,
    so these are excluded from the model step."""
    from alphagenome.data import genome
    try:
        client.query_variant(genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT),
                             requested_scorers=['AVI_SCORE'])
        return 'no scores'
    except ValueError as err:
        return 'reference mismatch' if 'reference base' in str(err) else 'not in Atlas'
    except Exception:  # noqa: BLE001
        return 'query failed'


def run_atlas(args, variants, scorers, targets, cache):
    from alphagenome.atlas import atlas
    client = atlas.create(load_api_key(args.api_key_file), timeout=180)
    size, workers = args.batch_size or 100, args.max_workers or 10
    snv = variants[(variants.REF.str.len() == 1) & (variants.ALT.str.len() == 1)]
    frames = []
    for gene, group in snv.groupby('gene', sort=True):
        for i, batch in chunks(group, size):
            ids = ','.join(batch.variant_id)
            path = cache / f"{gene}_{hashlib.md5(ids.encode()).hexdigest()[:12]}.parquet"
            if not path.exists():
                atlas_batch(client, batch, gene, scorers, targets, workers).to_parquet(path, index=False)
            frames.append(pd.read_parquet(path))
        print(f"{gene}: {len(group):,} SNVs", file=sys.stderr, flush=True)
    out = pd.concat(frames, ignore_index=True)
    scored = out.filter(like='|q').notna().any(axis=1)
    missing = variants[~variants.variant_id.isin(out.variant_id[scored])].copy()
    missing['atlas_note'] = [atlas_reason(client, r) for r in missing.itertuples()]
    if args.missing:
        missing.to_csv(args.missing, sep='\t', index=False)
    return out[scored], missing


# --------------------------------------------------------------------------------------
# score: model API
# --------------------------------------------------------------------------------------
def run_model(args, variants, scorers, targets, cache):
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    model = dna_client.create(load_api_key(args.api_key_file))
    recommended = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected = [recommended[s] for s in scorers]
    names = {recommended[s]: s for s in scorers}
    length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{args.sequence_length}']
    size, workers = args.batch_size or 10, args.max_workers or 5
    frames, n_batches = [], -(-len(variants) // size)
    for i, batch in chunks(variants, size):
        ids = ','.join(batch.variant_id)
        path = cache / f"model_{hashlib.md5(ids.encode()).hexdigest()[:12]}.parquet"
        if not path.exists():
            vs = [genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT, name=r.variant_id)
                  for r in batch.itertuples()]
            results = with_retries(
                lambda: model.score_variants([v.reference_interval.resize(length) for v in vs], vs,
                                             variant_scorers=selected, progress_bar=False,
                                             max_workers=workers),
                f'model batch {i + 1}/{n_batches}')
            per_scorer = {}
            for v, gene, per_variant in zip(vs, batch.gene, results):
                key = str(genome.Variant(v.chromosome, v.position, v.reference_bases, v.alternative_bases))
                for ad in per_variant:
                    scorer = names[ad.uns['variant_scorer']]
                    per_scorer.setdefault(scorer, []).append(
                        reduce_scorer(scorer, ad, [key] * ad.n_obs, [gene] * ad.n_obs, targets))
            frame = [pd.concat(f) for f in per_scorer.values()]
            summarise(frame, batch, 'model').to_parquet(path, index=False)
        frames.append(pd.read_parquet(path))
        print(f"model batch {i + 1}/{n_batches}", file=sys.stderr, flush=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['variant_id'])

def run_score(args):
    variants = read_variants([args.input])
    default = ATLAS_SCORERS if args.backend == 'atlas' else MODEL_SCORERS
    scorers = [s for s in (args.scorers or ','.join(default)).split(',') if s]
    if args.backend == 'model':
        scorers = [s for s in scorers if s != 'AVI_SCORE']
    targets = [t for t in args.target_ontology.split(',') if t]
    tag = hashlib.md5(f"{args.backend}{scorers}{targets}{args.sequence_length}".encode()).hexdigest()[:8]
    cache = (args.cache or args.output.parent / f'.{args.output.stem}_cache') / tag
    cache.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if variants.empty:
        out = pd.DataFrame(columns=['variant_id', 'gene', 'CHROM', 'POS', 'REF', 'ALT', 'backend'])
        if args.missing:
            variants.to_csv(args.missing, sep='\t', index=False)
    elif args.backend == 'atlas':
        out, missing = run_atlas(args, variants, scorers, targets, cache)
        print(f"Atlas: {len(out):,}/{len(variants):,} SNVs scored, {len(missing):,} missing",
              file=sys.stderr)
    else:
        out = run_model(args, variants, scorers, targets, cache)
        print(f"Model: {len(out):,}variants scored", file=sys.stderr)
    out.to_parquet(args.output, index=False)


# --------------------------------------------------------------------------------------
# merge
# --------------------------------------------------------------------------------------
def run_merge(args):
    frames = [pd.read_parquet(p) for p in args.summaries]
    df = pd.concat([f for f in frames if len(f)], ignore_index=True)
    df = df.drop_duplicates('variant_id', keep='first')   # Atlas first, model for the rest
    ann = pd.read_csv(args.annotated, dtype={'chrom': str}).drop_duplicates('Variant_ID')
    cols = [c for c in ['VariantID', 'type', 'region', 'detail', 'Consequence', 'Protein_position',
                        'Amino_acids', 'exon_distance'] if c in ann]
    df = df.merge(ann[cols].rename(columns={'VariantID': 'variant_id'}), on='variant_id', how='left')
    lead = ['variant_id', 'gene', 'CHROM', 'POS', 'REF', 'ALT', 'type', 'region', 'detail', 'backend']
    lead = [c for c in lead if c in df]
    df = df[lead + [c for c in df.columns if c not in lead]].sort_values(['gene', 'POS'])
    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out_parquet, index=False)
    df.to_csv(args.out_tsv, sep='\t', index=False)
    print(f"{len(df):,} variants with AlphaGenome modality scores "
          f"({(df.backend == 'atlas').sum():,} Atlas, {(df.backend == 'model').sum():,} model)",
          file=sys.stderr)


def main(argv=None):
    args = parse_args(argv)
    {'select': run_select, 'score': run_score, 'merge': run_merge}[args.cmd](args)


if __name__ == '__main__':
    main()