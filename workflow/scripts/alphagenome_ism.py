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

CONTEXT_1MB = 2**20  # dna_client.SEQUENCE_LENGTH_1MB

# Track-level Atlas scorers.
ATLAS_SCORERS = ['AVI_SCORE', 'RNA_SEQ', 'SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS',
                 'POLYADENYLATION', 'CAGE', 'PROCAP', 'DNASE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF',
                 'CONTACT_MAPS']

sys.path.insert(0, str(Path(__file__).resolve().parent))




def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    sub = p.add_subparsers(dest='cmd', required=True)

    s = sub.add_parser('scan')
    s.add_argument('GENE', help='Gene whose TSS is in center of query')
    s.add_argument('outdir')
    s.add_argument('--gtf')
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


def gene_context(gene: str, gtf_path: str, context: int = CONTEXT_1MB) -> genome.Interval:
    """1 Mb interval centred on the TSS of `gene`. Swap in your own annotation."""
    cols = ["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"]
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", names=cols,
                      dtype={"chrom": str}, low_memory=False)
    hit = gtf[(gtf.feat == "gene") & gtf.attr.str.contains(f'gene_name "{gene}:')]
    if hit.empty:
        raise ValueError(f"{gene} not found in {gtf_path}")
    g = hit.iloc[0]
    tss = int(g.start) if g.strand == "+" else int(g.end)  # GTF is 1-based
    half = context // 2
    # genome.Interval is 0-based half-open.
    return genome.Interval(g.chrom, max(0, tss - half), max(0, tss - half) + context)





def main(argv=None):
    args = parse_args(argv)
    {'scan': run_scan, 'plot': run_plot}[args.cmd](args)

if __name__ == '__main__':
    main()