"""
alphagenome_ism.py
================================

run ISM around gene of interest from the alphagenome atlas using the alphagenome atlas API
"""

from __future__ import annotations

import sys
import numpy as np
import pandas as pd
import argparse
from pathlib import Path

from alphagenome.atlas import atlas
from alphagenome.data import genome
from alphagenome_utils import load_api_key

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vep_report import load_gtf

CONTEXT_1MB = 2**20  # dna_client.SEQUENCE_LENGTH_1MB

# Track-level Atlas scorers.
ATLAS_SCORERS = ['AVI_SCORE', 'RNA_SEQ', 'SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS',
                 'POLYADENYLATION', 'CAGE', 'PROCAP', 'DNASE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF',
                 'CONTACT_MAPS']


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('scan')
    s.add_argument('gene', help='Gene whose TSS is in center of query')
    s.add_argument('outdir')
    s.add_argument('--gtf', required=True, help="GENCODE hg38 GTF as feather (path or URL)")
    s.add_argument('--window', type=int, default=10_000, help='window to query')
    s.add_argument('--target-ontology', default='CL:0002327,UBERON:0008367')
    s.add_argument('--tile-size', type=int, default=2048)
    s.add_argument('--max-workers', type=int, default=10)
    s.add_argument('--scorers', default=None)
    s.add_argument('--api-key-file', default='resources/API_KEY')

    s = sub.add_parser('plot')
    s.add_argument('outdir')
    s.add_argument('figdir')
    s.add_argument('--gtf')
    return p.parse_args(argv)


def gene_context(gene: str, gtf_path: str, width: int = CONTEXT_1MB) -> dict:
    """
    Interval of width centred on the TSS of `gene`. Returns gene, gene_id, transcript_id,
    chrom, strand, tss, start and end.
    tss is 1-based; start/end are 0-based half-open (genome.Interval)
    """
    _, mane = load_gtf(gtf_path, [gene])    # mane: one row per gene, indexed by gene_name
    t = mane.loc[gene]                      # has Chromosome, Strand, tss (1-based), gene_id, transcript_id, Start, End
    tss = int(t.tss)
    half = width // 2
    start = max(0, (tss-1)-half)
    context_dict = {'gene': gene, 'gene_id': t.gene_id, 'transcript_id': t.transcript_id,
                    'chrom': t.Chromosome, 'strand': t.Strand, 'tss': tss,
                    'start': start, 'end': start + width}
    return context_dict


def scorer_plan(metadata, scorers, targets) -> dict:
    """
    Decide, before querying, what each scorer will return.

    metadata: client.scorer_metadata(); scorers: scorer names; targets: ontology CURIEs
    ([] = all tracks).

    Returns {scorer: {'filter': bool, 'tracks': DataFrame, 'signed': bool}}:
      filter  - query this scorer with ontology_terms=targets (False = all tracks)
      tracks  - the tracks the Atlas will return, in column order, numbered 0..n-1
      signed  - scores are ALT - REF with a sign (False = magnitude only)

    A scorer with no tracks for `targets` (e.g. PROCAP, CONTACT_MAPS for breast) falls
    back to all its tracks, instead of being queried and returning nothing.

    Usage (e.g. in ism_dev.ipynb):
        client = atlas.create(load_api_key('resources/API_KEY'))
        plan = scorer_plan(client.scorer_metadata(), ATLAS_SCORERS,
                           ['CL:0002327', 'UBERON:0008367'])
        {s: (len(p['tracks']), p['filter'], p['signed']) for s, p in plan.items()}
        # expect e.g. DNASE (2, True, True), PROCAP (12, False, True),
        #             CONTACT_MAPS (28, False, False)
    """
    plan = {}
    for s in scorers:
        # 1. Unknown scorer name: fail now, not an hour into the scan
        if s not in metadata:
            raise ValueError(f"Unknown Atlas scorer {s!r}; available: {sorted(metadata)}")
        tracks = metadata[s].track_metadata.reset_index(drop=True)

        # 2. Filtering only applies if targets were given and the scorer has tissue tracks
        #    (AVI_SCORE and SPLICE_SITES have no ontology_curie column)
        use_filter = bool(targets) and 'ontology_curie' in tracks.columns

        # 3. Filter; fall back to all tracks if none match
        if use_filter:
            hit = tracks[tracks['ontology_curie'].isin(targets)].reset_index(drop=True)
            if hit.empty:
                print(f"{s}: no tracks for {','.join(targets)}; using all {len(tracks)} tracks",
                      file=sys.stderr)
                use_filter = False
            else:
                tracks = hit

        # 4. Store
        plan[s] = {'filter': use_filter, 'tracks': tracks, 'signed': metadata[s].is_signed}
    return plan




def run_scan(args):
    """PLACEHOLDER, for running the actual scan"""
    targets = [t for t in args.target_ontology.split(',') if t]
    w = gene_context(args.gene, args.gtf, args.window)

    client = atlas.create(load_api_key(KEY))
    plan = scorer_plan(client.scorer_metadata(), ATLAS_SCORERS, targets)



    iv = genome.Interval(w)
    ads = {}
    for f in (True, False):
        group = [s for s, p in plan.items() if p['filter'] == f]
        ads.update(client.query_interwal(iv, request_scorers=group,
                                         ontology_terms=targets if f else None,
                                         gene_names=['BRCA1'], progress_bar=False))
        






    print(args)


def run_plot(args):
    """PLACEHOLDER, for plotting the genome track"""
    print(args)







def main(argv=None):
    args = parse_args(argv)
    {'scan': run_scan, 'plot': run_plot}[args.cmd](args)

if __name__ == '__main__':
    main()