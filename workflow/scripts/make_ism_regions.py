#!/usr/bin/env python3
"""
make_ism_regions.py
===================

One ISM window per target gene: the MANE Select transcript extended by
``--flank`` bp on both sides (or, with ``--around-tss``, the MANE TSS +/- flank),
clipped at 0. Output (0-based, half-open): ``chrom, start, end, gene, strand, tss``.

Usage
-----
::

    python make_ism_regions.py GTF_FEATHER OUT_TSV --genes G1,G2 [--flank 100000] [--around-tss]
"""

import argparse

import numpy as np
import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('gtf', help='GENCODE hg38 GTF as feather (path or URL), as in vep_report.py')
    p.add_argument('out_tsv')
    p.add_argument('--genes', required=True)
    p.add_argument('--flank', type=int, default=100_000)
    p.add_argument('--around-tss', action='store_true')
    args = p.parse_args(argv)

    genes = args.genes.split(',')
    gtf = pd.read_feather(args.gtf, columns=['Chromosome', 'Start', 'End', 'Strand', 'Feature',
                                             'gene_name', 'tag'])
    mane = gtf[(gtf.Feature == 'transcript') & gtf.gene_name.isin(genes)
               & gtf.tag.str.contains('MANE_Select', na=False)].drop_duplicates('gene_name')
    missing = sorted(set(genes) - set(mane.gene_name))
    if missing:
        raise ValueError(f"No MANE Select transcript for: {missing}")
    tss = np.where(mane.Strand == '-', mane.End, mane.Start + 1)
    lo, hi = (tss, tss) if args.around_tss else (mane.Start.values, mane.End.values)
    out = pd.DataFrame({
        'chrom': mane.Chromosome.values,
        'start': np.maximum(lo - args.flank, 0),
        'end': hi + args.flank,
        'gene': mane.gene_name.values,
        'strand': mane.Strand.values,
        'tss': tss,
    }).sort_values('gene')
    out.to_csv(args.out_tsv, sep='\t', index=False)
    print(f"{len(out)} regions, {(out.end - out.start).sum() / 1e6:.1f} Mb in total "
          f"(~{3 * (out.end - out.start).sum() / 1e6:.0f} M SNVs)")


if __name__ == '__main__':
    main()
