#!/usr/bin/env python3
"""
query_alphagenome_atlas.py
==========================

Fetch precomputed AlphaGenome scores from the AlphaGenome Atlas (``alphagenome.atlas``),
as in ``notebooks/alphagenome_atlas_BRCA1.ipynb``, for a list of SNVs or for every
SNV in a set of regions (in silico mutagenesis, ISM).

The Atlas returns stored scores, so it is far cheaper than running the model:
use it for SNVs and the model API (score_variants_alphagenome.py) for indels and
for anything the Atlas does not hold.

Modes
-----
``variants VARIANTS_TSV OUT``
    VARIANTS_TSV has ``variant_id, CHROM, POS, REF, ALT, gene`` (from
    select_alphagenome_variants.py). Variants are queried per target gene, so
    gene-level scorers (RNA_SEQ, SPLICE_SITES, ...) report the effect on that gene.
    Variants the Atlas has no scores for are written to ``--missing``.
``interval REGIONS_TSV OUT``
    REGIONS_TSV has ``chrom, start, end, gene`` (0-based, half-open); every SNV in
    each region is fetched in ``--chunk`` bp pieces, effects reported for ``gene``.

Output
------
``OUT`` (parquet, or TSV if it ends in .tsv/.tsv.gz): one row per variant (and
gene), columns ``variant_id, gene, chrom, pos, ref, alt`` then one column per
scorer track, named ``<SCORER>:<track name>[:<strand>]``, plus ``<column>|q`` for
the calibrated quantile where the Atlas provides one.

Pieces are cached under ``--cache`` (one parquet per gene / chunk), so a rerun
after an interruption only fetches what is missing.

Usage
-----
::

    python query_alphagenome_atlas.py variants snvs.tsv atlas.parquet --missing missing.tsv
        [--scorers AVI_SCORE,RNA_SEQ,...] [--ontology-terms CL:0002327,...] [--cache DIR]
    python query_alphagenome_atlas.py interval regions.tsv ism.parquet --chunk 4000
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alphagenome_utils import chunks, load_api_key, read_variants, with_retries  # noqa: E402

DEFAULT_SCORERS = 'AVI_SCORE,RNA_SEQ,SPLICE_SITES,DNASE,CAGE,CHIP_HISTONE'
DEFAULT_ONTOLOGY = 'CL:0002327'   # mammary epithelial cell


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('mode', choices=['variants', 'interval'])
    p.add_argument('input')
    p.add_argument('output', type=Path)
    p.add_argument('--missing', type=Path, default=None,
                   help='variants mode: write variants without Atlas scores here')
    p.add_argument('--scorers', default=DEFAULT_SCORERS)
    p.add_argument('--ontology-terms', default=DEFAULT_ONTOLOGY,
                   help="Comma-separated ontology CURIEs; '' for all biosamples")
    p.add_argument('--cache', type=Path, default=None)
    p.add_argument('--batch-size', type=int, default=200, help='variants mode: variants per query')
    p.add_argument('--chunk', type=int, default=4000, help='interval mode: bp per query')
    p.add_argument('--gene', default=None, help='interval mode: only the region(s) of this gene')
    p.add_argument('--max-workers', type=int, default=10)
    p.add_argument('--api-key-file', default=None)
    p.add_argument('--timeout', type=float, default=120)
    return p.parse_args(argv)


def tidy_atlas(result, gene):
    """One row per variant, one column per scorer track (raw score, and ``|q`` quantile)."""
    frames = []
    for scorer, ad in result.items():
        if ad.n_obs == 0 or ad.obs is None or 'variant' not in ad.obs:
            continue  # e.g. splice scores outside the gene body
        names = ad.var['name'].astype(str).str.replace('Histone ChIP-seq', '').str.strip()
        strands = (ad.var['strand'].map(lambda s: '' if s == '.' else f':{s}')
                   if 'strand' in ad.var else '')
        cols = [f'{scorer}:{n}{s}' for n, s in zip(names, strands)]
        df = pd.DataFrame(ad.X, columns=cols)
        if 'quantiles' in ad.layers:
            df = df.join(pd.DataFrame(ad.layers['quantiles'], columns=[c + '|q' for c in cols]))
        df['variant'] = ad.obs['variant'].astype(str).values
        if 'gene_name' in ad.obs:
            df = df[(ad.obs['gene_name'] == gene).values]
        df = df.loc[:, ~df.columns.duplicated()]
        frames.append(df.drop_duplicates('variant').set_index('variant'))
    if not frames:
        return pd.DataFrame(columns=['variant'])
    return pd.concat(frames, axis=1).reset_index()


def split_variant_key(df):
    parts = df['variant'].str.extract(
        r'^(?P<chrom>[^:]+):(?P<pos>\d+):(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$')
    return pd.concat([parts.astype({'pos': 'Int64'}), df.drop(columns='variant')], axis=1)


def write_table(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(('.tsv', '.tsv.gz')):
        df.to_csv(path, sep='\t', index=False)
    else:
        df.to_parquet(path, index=False)


class Atlas:
    def __init__(self, args):
        from alphagenome.atlas import atlas
        self.client = atlas.create(load_api_key(args.api_key_file), timeout=args.timeout)
        self.scorers = [s for s in args.scorers.split(',') if s]
        self.ontology = [t for t in args.ontology_terms.split(',') if t] or None
        self.max_workers = args.max_workers
        tag = hashlib.md5(f"{self.scorers}{self.ontology}".encode()).hexdigest()[:8]
        self.cache = (args.cache or args.output.parent / f'.{args.output.stem}_cache') / tag
        self.cache.mkdir(parents=True, exist_ok=True)

    def query_variants(self, variants, gene):
        return self.client.query_variants(
            variants, requested_scorers=self.scorers, ontology_terms=self.ontology,
            gene_names=[gene], progress_bar=False, max_workers=self.max_workers)

    def query_interval(self, interval, gene):
        return self.client.query_interval(
            interval, requested_scorers=self.scorers, ontology_terms=self.ontology,
            gene_names=[gene], progress_bar=False, max_workers=self.max_workers)


def run_variants(args, atlas):
    from alphagenome.data import genome
    df = read_variants([args.input])
    if 'gene' not in df:
        raise ValueError(f"{args.input} needs a 'gene' column (target gene per variant)")
    df = df[(df.REF.str.len() == 1) & (df.ALT.str.len() == 1)]
    frames, missing = [], []
    for gene, group in df.groupby('gene', sort=True):
        for i, batch in chunks(group, args.batch_size):
            ids = ','.join(batch.variant_id)
            path = atlas.cache / f"{gene}_{hashlib.md5(ids.encode()).hexdigest()[:12]}.parquet"
            if not path.exists():
                variants = [genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT, name=r.variant_id)
                            for r in batch.itertuples()]
                try:
                    result = with_retries(lambda: atlas.query_variants(variants, gene),
                                          f'{gene} batch {i}', attempts=3)
                    tidy = tidy_atlas(result, gene)
                except ValueError:
                    # A variant the Atlas does not hold fails the whole batch: retry one by one.
                    parts = []
                    for v in variants:
                        try:
                            parts.append(tidy_atlas(atlas.query_variants([v], gene), gene))
                        except ValueError as err:
                            print(f"not in Atlas: {v} ({str(err)[:120]})", file=sys.stderr)
                    tidy = pd.concat(parts, ignore_index=True) if parts else tidy_atlas({}, gene)
                tidy.to_parquet(path, index=False)
            tidy = pd.read_parquet(path)
            keys = {str(genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT)): r.variant_id
                    for r in batch.itertuples()}
            tidy['variant_id'] = tidy['variant'].map(keys)
            tidy['gene'] = gene
            frames.append(tidy)
            found = set(tidy.variant_id.dropna())
            missing.append(batch[~batch.variant_id.isin(found)])
        print(f"{gene}: {len(group)} SNVs queried", file=sys.stderr)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['variant'])
    out = split_variant_key(out)
    lead = ['variant_id', 'gene', 'chrom', 'pos', 'ref', 'alt']
    out = out[lead + [c for c in out.columns if c not in lead]]
    write_table(out, args.output)
    miss = pd.concat(missing, ignore_index=True) if missing else df.iloc[:0]
    if args.missing:
        miss.to_csv(args.missing, sep='\t', index=False)
    print(f"Atlas scores for {out.variant_id.nunique()}/{len(df)} SNVs; "
          f"{len(miss)} missing", file=sys.stderr)


def run_interval(args, atlas):
    from alphagenome.data import genome
    regions = pd.read_csv(args.input, sep='\t')
    if args.gene:
        regions = regions[regions.gene == args.gene]
        if regions.empty:
            raise ValueError(f"No region for gene {args.gene} in {args.input}")
    frames = []
    for r in regions.itertuples():
        for start in range(int(r.start), int(r.end), args.chunk):
            end = min(start + args.chunk, int(r.end))
            path = atlas.cache / f"{r.gene}_{r.chrom}_{start}_{end}.parquet"
            if not path.exists():
                interval = genome.Interval(r.chrom, start, end)
                result = with_retries(lambda: atlas.query_interval(interval, r.gene),
                                      f'{r.gene} {r.chrom}:{start}-{end}')
                tidy_atlas(result, r.gene).to_parquet(path, index=False)
            frames.append(pd.read_parquet(path).assign(gene=r.gene))
        print(f"{r.gene}: {r.chrom}:{r.start}-{r.end} done", file=sys.stderr)
    out = split_variant_key(pd.concat(frames, ignore_index=True))
    out['variant_id'] = (out.chrom + '_' + out.pos.astype(str) + '_' + out.ref + '_' + out.alt)
    lead = ['variant_id', 'gene', 'chrom', 'pos', 'ref', 'alt']
    out = out[lead + [c for c in out.columns if c not in lead]].sort_values(['gene', 'pos', 'alt'])
    write_table(out, args.output)


def main(argv=None):
    args = parse_args(argv)
    atlas = Atlas(args)
    (run_variants if args.mode == 'variants' else run_interval)(args, atlas)


if __name__ == '__main__':
    main()
