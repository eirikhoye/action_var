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


def run_scan(args):
    """PLACEHOLDER, for running the actual scan"""
    print(args)


def run_plot(args):
    """PLACEHOLDER, for plotting the genome track"""
    print(args)


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





def main(argv=None):
    args = parse_args(argv)
    {'scan': run_scan, 'plot': run_plot}[args.cmd](args)

if __name__ == '__main__':
    main()