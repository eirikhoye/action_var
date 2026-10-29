#!/usr/bin/env python3
"""
score_variants_alphagenome.py
=============================

Score variants with the AlphaGenome model API (``dna_client.score_variants``).

Use this for variants the AlphaGenome Atlas does not hold (indels, and SNVs the
Atlas has no scores for); SNVs are cheaper through query_alphagenome_atlas.py.

Input
-----
One or more TSVs with ``variant_id, CHROM, POS, REF, ALT`` and optionally
``gene`` (the target gene; empty files are fine). Duplicate IDs are dropped.

Method
------
* Variants are scored in batches of ``--batch-size`` with the recommended
  scorers selected by ``--scorers`` (or switched off with ``--no-<scorer>``),
  in a ``--sequence-length`` window centred on the variant.
* To keep the output manageable, each batch's tidy scores are filtered to the
  ``--ontology-terms`` biosamples (``''`` keeps all) and, for gene-level scorers,
  to the variant's ``gene`` (``--all-genes`` keeps every gene in the window).
* Each batch is cached as parquet under ``--cache``; a rerun only scores
  batches that are missing. Transient API errors are retried with back-off.

Output
------
Long table, one row per (variant, gene, track), as
``alphagenome.models.variant_scorers.tidy_scores`` (raw_score, quantile_score,
track and biosample metadata). Parquet unless the path ends in .csv/.csv.gz.

API key: ``--api-key-file``, else ``$ALPHAGENOME_API_KEY``, else
``resources/API_KEY``. Never commit the key.

Usage
-----
::

    python score_variants_alphagenome.py OUT INPUT_TSV [INPUT_TSV ...]
        [--scorers RNA_SEQ,SPLICE_SITES,...] [--ontology-terms CL:0002327,...]
        [--sequence-length 1MB] [--batch-size 20] [--cache DIR]
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alphagenome_utils import chunks, load_api_key, read_variants, with_retries  # noqa: E402

SCORER_KEYS = ['RNA_SEQ', 'CAGE', 'PROCAP', 'ATAC', 'DNASE', 'CHIP_HISTONE', 'CHIP_TF',
               'POLYADENYLATION', 'SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS']
DEFAULT_ONTOLOGY = 'CL:0002327'   # mammary epithelial cell


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('output', type=Path)
    p.add_argument('inputs', nargs='+', help='variant TSV(s)')
    p.add_argument('--organism', default='human', choices=['human', 'mouse'])
    p.add_argument('--sequence-length', '--sequence_length', dest='sequence_length',
                   default='1MB', choices=['16KB', '100KB', '500KB', '1MB'])
    p.add_argument('--scorers', default=','.join(SCORER_KEYS),
                   help='Comma-separated recommended scorers (default: %(default)s)')
    for key in SCORER_KEYS:   # backwards compatible --no-<scorer> switches
        p.add_argument(f"--no-{key.lower().replace('_', '-')}", dest='disabled',
                       action='append_const', const=key, help=f'Disable the {key} scorer')
    p.add_argument('--ontology-terms', default=DEFAULT_ONTOLOGY,
                   help="Comma-separated ontology CURIEs to keep; '' keeps all")
    p.add_argument('--all-genes', action='store_true',
                   help="Keep gene-level scores for every gene, not only the variant's gene")
    p.add_argument('--batch-size', type=int, default=20)
    p.add_argument('--max-workers', type=int, default=5)
    p.add_argument('--cache', type=Path, default=None)
    p.add_argument('--api-key-file', default=None)
    return p.parse_args(argv)


def select_scorers(args, organism):
    from alphagenome.models import dna_client, variant_scorers
    wanted = [s for s in args.scorers.split(',') if s and s not in (args.disabled or [])]
    unknown = set(wanted) - set(variant_scorers.RECOMMENDED_VARIANT_SCORERS)
    if unknown:
        raise ValueError(f"Unknown scorer(s) {sorted(unknown)}; "
                         f"available: {list(variant_scorers.RECOMMENDED_VARIANT_SCORERS)}")
    selected = []
    for key in wanted:
        scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[key]
        supported = organism.value in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        if scorer.requested_output == dna_client.OutputType.PROCAP and \
                organism == dna_client.Organism.MUS_MUSCULUS:
            supported = False
        if supported:
            selected.append(scorer)
        else:
            print(f"Excluding {key}: not supported for {organism}", file=sys.stderr)
    return selected


def filter_scores(df, batch, ontology, all_genes):
    """Keep the requested biosamples, and gene-level rows for the variant's own gene."""
    if df is None or df.empty:
        return pd.DataFrame()
    df = df.copy()
    df['variant_id'] = df['variant_id'].astype(str)
    if ontology and 'ontology_curie' in df:
        df = df[df.ontology_curie.isin(ontology)]
    if not all_genes and 'gene' in batch and 'gene_name' in df:
        names = {str(v): g for v, g in zip(batch.variant_key, batch.gene)}
        target = df.variant_id.map(names)
        df = df[df.gene_name.isna() | (df.gene_name == target) | target.isna()]
    for c in df.columns:   # parquet needs plain types (Interval / Variant objects -> str)
        if df[c].dtype == object:
            df[c] = df[c].where(df[c].isna(), df[c].astype(str))
    return df


def main(argv=None):
    args = parse_args(argv)
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    variants = read_variants(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if variants.empty:
        print("No variants to score", file=sys.stderr)
        pd.DataFrame(columns=['variant_id']).to_parquet(args.output, index=False)
        return

    organism = {'human': dna_client.Organism.HOMO_SAPIENS,
                'mouse': dna_client.Organism.MUS_MUSCULUS}[args.organism]
    scorers = select_scorers(args, organism)
    length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{args.sequence_length}']
    ontology = [t for t in args.ontology_terms.split(',') if t]
    model = dna_client.create(load_api_key(args.api_key_file))

    tag = hashlib.md5(f"{args.scorers}{args.disabled}{length}{ontology}{args.all_genes}"
                      .encode()).hexdigest()[:8]
    cache = (args.cache or args.output.parent / f'.{args.output.stem}_cache') / tag
    cache.mkdir(parents=True, exist_ok=True)

    frames = []
    n_batches = -(-len(variants) // args.batch_size)
    for i, batch in chunks(variants, args.batch_size):
        ids = ','.join(batch.variant_id)
        path = cache / f"batch_{hashlib.md5(ids.encode()).hexdigest()[:12]}.parquet"
        if not path.exists():
            vs = [genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT, name=r.variant_id)
                  for r in batch.itertuples()]
            intervals = [v.reference_interval.resize(length) for v in vs]
            batch = batch.assign(variant_key=[str(v) for v in vs])
            scores = with_retries(
                lambda: model.score_variants(intervals, vs, variant_scorers=scorers,
                                             organism=organism, progress_bar=False,
                                             max_workers=args.max_workers),
                f'batch {i + 1}/{n_batches}')
            tidy = variant_scorers.tidy_scores(scores)
            filter_scores(tidy, batch, ontology, args.all_genes).to_parquet(path, index=False)
        frames.append(pd.read_parquet(path))
        print(f"batch {i + 1}/{n_batches} done", file=sys.stderr)

    out = pd.concat(frames, ignore_index=True)
    # tidy_scores names variants "chr:pos:ref>alt"; add the input IDs back
    keys = {str(genome.Variant(r.CHROM, int(r.POS), r.REF, r.ALT)): r.variant_id
            for r in variants.itertuples()}
    out = out.rename(columns={'variant_id': 'variant'})
    out.insert(0, 'variant_id', out.variant.map(keys))
    if str(args.output).endswith(('.csv', '.csv.gz')):
        out.to_csv(args.output, index=False)
    else:
        out.to_parquet(args.output, index=False)
    print(f"{out.variant_id.nunique()}/{len(variants)} variants scored, {len(out)} rows",
          file=sys.stderr)


if __name__ == '__main__':
    main()
