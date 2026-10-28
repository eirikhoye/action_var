#!/usr/bin/env python3
"""
vep_report.py
=============

Summarise a VEP-annotated VCF per target gene and region, as an HTML report.

Purpose
-------
Script version of ``workflow/notebooks/analyze_new_bridges34_variants.ipynb``,
run by the ``vep_report`` rule after ``annotate_vep``. It answers: which target
genes do the variants fall in, and where in those genes (coding, splice, UTR,
intron, promoter, ...)?

Input
-----
``VEP_VCF``
    Output of annotate_vep.py: GRCh38, ``ID`` = VariantID, and a ``CSQ`` field
    that includes ``MANE_SELECT`` (VEP run with ``--mane``).

Method
------
1. Variants are read from the VCF records; duplicate IDs are dropped.
2. Each variant gets the consequence on the MANE Select transcript of a target
   gene. VEP reports genes within 5 kb; when two target genes qualify, the one
   the variant lies inside wins, then the nearest.
3. Region class, when VEP gives several terms: canonical splice site > coding
   > splice region > UTR > intron > upstream / downstream.
4. Variants further than 5 kb from every target gene are "distal" and assigned
   to the nearest target-gene TSS on the same chromosome. Variants on a
   chromosome with no target gene are dropped (and counted in the report).
5. Gene models come from the GENCODE hg38 GTF given with ``--gtf``
   (MANE Select transcripts).

Output
------
``REPORT_HTML``
    Self-contained HTML report (figures embedded).
``--annotated``
    Per-variant table: hg38 position, gene, region, consequence, distances.
``--figdir``
    The report's figures as PNG, plus ``tracks/<gene>_variant_track.png``.

Usage
-----
::

    python vep_report.py VEP_VCF REPORT_HTML --annotated CSV --figdir DIR
        --genes G1,G2,... --gtf PATH_OR_URL

Dependencies
------------
Python >= 3.9, pandas, numpy, matplotlib, pyarrow (+ fsspec for a URL GTF).

Author: Eirik Høye
"""

import argparse
import base64
import html
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Region classes in 5'->3' order, each with a fixed categorical colour.
REGIONS = {
    'upstream (1-5 kb)': '#2a78d6',
    'promoter (≤1 kb)': '#eb6834',
    "5' UTR": '#1baf7a',
    'coding': '#eda100',
    'splice site / region': '#e87ba4',
    'intron': '#008300',
    "3' UTR": '#6250d6',
    'downstream (≤5 kb)': '#e34948',
    'distal (>5 kb)': '#b5b4ad',
}

CODING_TERMS = {
    'missense_variant': 'missense', 'synonymous_variant': 'synonymous',
    'stop_gained': 'stop gained', 'frameshift_variant': 'frameshift',
    'start_lost': 'start / stop lost', 'stop_lost': 'start / stop lost',
    'inframe_insertion': 'in-frame indel', 'inframe_deletion': 'in-frame indel',
    'protein_altering_variant': 'other coding', 'stop_retained_variant': 'synonymous',
    'start_retained_variant': 'synonymous', 'coding_sequence_variant': 'other coding',
    'incomplete_terminal_codon_variant': 'other coding',
}
CODING_ORDER = ['stop gained', 'frameshift', 'start / stop lost', 'in-frame indel', 'missense',
                'other coding', 'synonymous']
SPLICE_REGION_TERMS = {'splice_region_variant', 'splice_donor_region_variant',
                       'splice_donor_5th_base_variant', 'splice_polypyrimidine_tract_variant'}
DETAIL_ORDER = CODING_ORDER + ['canonical splice site', 'splice region']

INK, SECONDARY, MUTED, GRID, AXIS, BLUE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#2a78d6'
plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.edgecolor': AXIS, 'axes.labelcolor': SECONDARY,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'xtick.labelcolor': SECONDARY, 'ytick.labelcolor': SECONDARY,
    'axes.spines.top': False, 'axes.spines.right': False,
    'font.size': 12, 'axes.titlesize': 13, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
    'figure.titlesize': 16, 'figure.titleweight': 'bold',
})


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="HTML report of VEP-annotated variants per gene and region",
        epilog="See the module docstring for the method.",
    )
    parser.add_argument("vep_vcf", type=Path, help="VEP output VCF (run with --mane)")
    parser.add_argument("report_html", type=Path, help="Output HTML report")
    parser.add_argument("--annotated", type=Path, required=True,
                        help="Output per-variant annotation table (CSV)")
    parser.add_argument("--figdir", type=Path, required=True, help="Output directory for PNG figures")
    parser.add_argument("--genes", required=True,
                        help="Comma-separated target genes (HGNC symbols, each with a MANE Select "
                             "transcript in the GTF)")
    parser.add_argument("--gtf", required=True,
                        help="GENCODE hg38 GTF as feather (path or URL)")
    return parser.parse_args(argv)


# --------------------------------------------------------------------------------------
# Reading and annotating
# --------------------------------------------------------------------------------------
def read_vep_vcf(path, genes):
    """Return (variants, csq): one row per VCF record, and the CSQ entries on
    MANE Select transcripts of ``genes``."""
    fields, records, csq_rows = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ'):
                fields = line.split('Format: ')[1].split('"')[0].split('|')
                if 'MANE_SELECT' not in fields:
                    raise ValueError(f"{path} has no MANE_SELECT in CSQ; rerun VEP with --mane")
                i_symbol, i_mane = fields.index('SYMBOL'), fields.index('MANE_SELECT')
                continue
            if line.startswith('#'):
                continue
            if fields is None:
                raise ValueError(f"{path} has no CSQ header line; is it a VEP VCF?")
            chrom, pos, vid, ref, alt, _, _, info = line.rstrip('\n').split('\t')[:8]
            records.append((vid, chrom, int(pos), ref, alt))
            csq = next((x[4:] for x in info.split(';') if x.startswith('CSQ=')), '')
            for entry in csq.split(','):
                values = entry.split('|')
                if len(values) == len(fields) and values[i_symbol] in genes and values[i_mane]:
                    csq_rows.append([vid] + values)
    variants = pd.DataFrame(records, columns=['VariantID', 'chrom', 'pos', 'Ref', 'Alt'])
    csq = pd.DataFrame(csq_rows, columns=['VariantID'] + fields)
    return variants, csq


def variant_type(variants):
    ref, alt = variants.Ref.str.len(), variants.Alt.str.len()
    return np.select([(ref == 1) & (alt == 1), ref > alt, ref < alt],
                     ['SNV', 'deletion', 'insertion'], 'MNV / complex')


def region_class(consequence, distance):
    terms = set(consequence.split('&'))
    if terms & {'splice_donor_variant', 'splice_acceptor_variant'}:
        return 'splice site / region', 'canonical splice site'
    coding = [CODING_TERMS[t] for t in terms if t in CODING_TERMS]
    if coding:
        return 'coding', min(coding, key=CODING_ORDER.index)
    if terms & SPLICE_REGION_TERMS and 'intron_variant' in terms:
        return 'splice site / region', 'splice region'
    if '5_prime_UTR_variant' in terms:
        return "5' UTR", None
    if '3_prime_UTR_variant' in terms:
        return "3' UTR", None
    if 'intron_variant' in terms:
        return 'intron', None
    if 'upstream_gene_variant' in terms:
        return ('promoter (≤1 kb)' if int(distance) <= 1000 else 'upstream (1-5 kb)'), None
    if 'downstream_gene_variant' in terms:
        return 'downstream (≤5 kb)', None
    return 'other', consequence


def load_gtf(path, genes):
    """GENCODE GTF (Start 0-based, End 1-based inclusive) and the MANE transcripts of ``genes``."""
    gtf = pd.read_feather(path)
    gtf = gtf[gtf.Feature.isin(['gene', 'transcript', 'exon', 'CDS'])]
    mane = (gtf[(gtf.Feature == 'transcript') & gtf.tag.str.contains('MANE_Select', na=False)
                & gtf.gene_name.isin(genes)].set_index('gene_name'))
    missing = sorted(set(genes) - set(mane.index))
    if missing:
        raise ValueError(f"No MANE Select transcript in {path} for: {missing}")
    mane['tss'] = np.where(mane.Strand == '-', mane.End, mane.Start + 1)
    return gtf, mane


def annotate(variants, csq, mane):
    """Assign each variant to one gene and region class; returns (ann, n_dropped)."""
    csq = csq.copy()
    csq['dist'] = pd.to_numeric(csq.DISTANCE.replace('', '0'))
    csq = csq.sort_values(['VariantID', 'dist']).drop_duplicates('VariantID')
    csq[['region', 'detail']] = [region_class(c, d) for c, d in zip(csq.Consequence, csq.dist)]

    ann = variants.merge(
        csq[['VariantID', 'SYMBOL', 'Feature', 'MANE_SELECT', 'Consequence', 'region', 'detail',
             'Protein_position', 'Amino_acids']],
        on='VariantID', how='left',
    ).rename(columns={'SYMBOL': 'gene', 'Feature': 'refseq_transcript'})

    # Distal variants: nearest target-gene TSS on the same chromosome.
    no_gene_on_chrom = ann.gene.isna() & ~ann.chrom.isin(mane.Chromosome)
    ann = ann[~no_gene_on_chrom].copy()
    for i, v in ann[ann.gene.isna()].iterrows():
        same = mane[mane.Chromosome == v.chrom]
        ann.loc[i, ['gene', 'region']] = (same.tss - v.pos).abs().idxmin(), 'distal (>5 kb)'

    # Signed distance to the TSS along the direction of transcription.
    t = mane.loc[ann.gene]
    ann['tss_distance'] = (ann.pos.values - t.tss.values) * np.where(t.Strand.values == '-', -1, 1)
    return ann, int(no_gene_on_chrom.sum())


def add_exon_distance(ann, gtf, mane):
    """Distance from each intronic / splice-region variant to the nearest MANE exon boundary."""
    exons = gtf[(gtf.Feature == 'exon') & gtf.transcript_id.isin(mane.transcript_id)]
    intronic = ann[ann.region.isin(['intron', 'splice site / region'])]
    ann['exon_distance'] = np.nan
    for gene, df in intronic.groupby('gene'):
        ex = exons[exons.gene_name == gene]
        starts, ends = ex.Start.values + 1, ex.End.values
        p = df.pos.values[:, None]
        ann.loc[df.index, 'exon_distance'] = np.minimum(np.abs(p - starts), np.abs(p - ends)).min(axis=1)


# --------------------------------------------------------------------------------------
# Tables and figures
# --------------------------------------------------------------------------------------
def region_tables(df, used):
    overall = df.region.value_counts().reindex(used).fillna(0).astype(int).to_frame('variants')
    overall['%'] = (100 * overall.variants / len(df)).round(1)
    by_type = pd.crosstab(df.region, df.type).reindex(used).fillna(0).astype(int)
    by_gene = pd.crosstab(df.gene, df.region).reindex(columns=used, fill_value=0)
    by_gene = by_gene.loc[by_gene.sum(axis=1).sort_values(ascending=False).index]
    return overall, by_type, by_gene


def detail_counts(df):
    return df[df.detail.isin(DETAIL_ORDER)].detail.value_counts().reindex(DETAIL_ORDER).fillna(0).astype(int)


def style_hbar(*axes):
    for a in axes:
        a.grid(axis='x', color=GRID, lw=0.6)
        a.set_axisbelow(True)
        a.tick_params(axis='y', length=0)


def fig_per_gene(by_gene, used, n_variants):
    fig, (ax, ax_pct) = plt.subplots(
        1, 2, figsize=(14, max(4, 0.29 * len(by_gene))), sharey=True,
        gridspec_kw={'width_ratios': [1.1, 1], 'wspace': 0.06})
    y = np.arange(len(by_gene))[::-1]
    for a, data in [(ax, by_gene), (ax_pct, by_gene.div(by_gene.sum(axis=1), axis=0) * 100)]:
        left = np.zeros(len(data))
        for region in used:
            a.barh(y, data[region], left=left, color=REGIONS[region], height=0.75, label=region,
                   edgecolor='white', lw=1)
            left += data[region].values
    style_hbar(ax, ax_pct)
    ax.set_yticks(y, by_gene.index, style='italic')
    ax.set_xlabel('variants')
    ax_pct.set_xlabel("% of the gene's variants")
    ax_pct.set_xlim(0, 100)
    ax.set_title('Number of variants')
    ax_pct.set_title('Share by region')
    ax.legend(ncols=4, frameon=False, loc='lower left', bbox_to_anchor=(0, 1.04), fontsize=11)
    fig.suptitle(f'{n_variants:,} variants across {len(by_gene)} genes, by region of the MANE transcript',
                 x=0.125, ha='left', y=0.99)
    return fig


def fig_region_and_consequence(overall, detail, title):
    overall = overall[overall.variants > 0]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={'width_ratios': [1.2, 1], 'wspace': 0.6})
    ax.barh(overall.index[::-1], overall.variants[::-1], color=[REGIONS[r] for r in overall.index[::-1]],
            height=0.7)
    ax.bar_label(ax.containers[0], labels=[f'{n:,}  ({p}%)' for n, p in
                                           zip(overall.variants[::-1], overall['%'][::-1])],
                 padding=4, fontsize=11, color=SECONDARY)
    ax.set_title(title)
    ax.set_xlabel('variants')
    ax.set_xlim(0, overall.variants.max() * 1.3)

    colors = [REGIONS['coding']] * len(CODING_ORDER) + [REGIONS['splice site / region']] * 2
    ax2.barh(detail.index[::-1], detail[::-1], color=colors[::-1], height=0.7)
    ax2.bar_label(ax2.containers[0], labels=[f'{n:,}' for n in detail[::-1]], padding=4, fontsize=11,
                  color=SECONDARY)
    ax2.set_title('Coding and splice variants by consequence')
    ax2.set_xlabel('variants')
    ax2.set_xlim(0, max(detail.max(), 1) * 1.2)
    style_hbar(ax, ax2)
    return fig


def fig_exon_distance(d, title):
    """Cumulative distribution of distances to the nearest exon."""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(d, np.arange(1, len(d) + 1) / len(d), color=BLUE, lw=2)
    for x in (50, 500, 5000):
        frac = (d <= x).mean()
        ax.annotate(f'{frac:.0%} within {x:,} bp', (x, frac), xytext=(8, -14), textcoords='offset points',
                    fontsize=11, color=INK)
        ax.plot(x, frac, 'o', color=INK, ms=5)
    ax.set_xscale('log')
    ax.set_xlabel('distance to nearest MANE exon (bp, log scale)')
    ax.set_ylabel('cumulative fraction')
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title)
    return fig


def draw_gene_models(ax, gtf, chrom, start, end, target):
    """MANE (or longest) transcript of each protein-coding gene in the window, packed into rows."""
    genes = gtf[(gtf.Feature == 'gene') & (gtf.Chromosome == chrom) & (gtf.End > start)
                & (gtf.Start < end) & (gtf.gene_type == 'protein_coding')
                & ~gtf.gene_name.str.startswith('ENSG')]
    genes = pd.concat([genes[genes.gene_name == target], genes[genes.gene_name != target]])
    row_end, min_width = [], (end - start) / 1500
    for g in genes.itertuples():
        tx = gtf[(gtf.gene_id == g.gene_id) & (gtf.Feature == 'transcript')]
        is_mane = tx.tag.str.contains('MANE_Select', na=False)
        tx = tx[is_mane].iloc[0] if is_mane.any() else tx.loc[(tx.End - tx.Start).idxmax()]
        lo, hi = max(tx.Start, start), min(tx.End, end)
        if lo >= hi:  # gene overlaps the window but its transcript does not
            continue
        row = 0 if g.gene_name == target else next(
            (i for i, e in enumerate(row_end) if i > 0 and lo > e), max(len(row_end), 1))
        row_end += [-np.inf] * (row + 1 - len(row_end))
        row_end[row] = hi + (end - start) * 0.12
        y, color = -row, INK if g.gene_name == target else MUTED
        ax.plot([lo, hi], [y, y], color=color, lw=1)
        parts = gtf[gtf.transcript_id == tx.transcript_id]
        for feature, h in [('exon', 0.15), ('CDS', 0.3)]:
            for f in parts[parts.Feature == feature].itertuples():
                a, b = max(f.Start, start), min(f.End, end)
                if a < b:  # at least ~1 px wide, so small exons of long genes stay visible
                    w = max(b - a, min_width)
                    ax.add_patch(plt.Rectangle(((a + b - w) / 2, y - h), w, 2 * h, color=color, lw=0))
        arrow = '←' if tx.Strand == '-' else '→'
        label = dict(fontsize=11, color=color, style='italic', clip_on=True)
        if g.gene_name == target:  # above its own row, which nothing else shares
            ax.text(lo, y + 0.4, f'{g.gene_name} {arrow}', ha='left', va='bottom', weight='bold', **label)
        else:  # to the right of the gene, inside the space reserved by the packing
            ax.text(hi, y, f' {g.gene_name} {arrow}', ha='left', va='center', **label)
    ax.set_ylim(-max(len(row_end), 1) + 0.5, 1.2)
    ax.set_yticks([])
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    ax.tick_params(bottom=False)


def fig_gene_track(ann, gtf, mane, gene, n_bins=200):
    df = ann[ann.gene == gene]
    t = mane.loc[gene]
    lo, hi = min(df.pos.min(), t.Start), max(df.pos.max(), t.End)
    pad = int((hi - lo) * 0.03) + 200
    start, end = lo - pad, hi + pad
    present = [r for r in REGIONS if (df.region == r).any()]

    fig, (ax_d, ax_v, ax_g) = plt.subplots(
        3, 1, figsize=(14, 3.2 + 0.32 * len(present)), sharex=True,
        gridspec_kw={'height_ratios': [1.2, 0.32 * len(present), 1.1], 'hspace': 0.15},
    )
    bin_bp = max(int(np.ceil((end - start) / n_bins / 50) * 50), 50)
    edges = np.arange(start, end + bin_bp, bin_bp)
    ax_d.hist(df.pos, bins=edges, color=BLUE, lw=0)
    ax_d.set_ylabel(f'variants\nper {bin_bp:,} bp', fontsize=11)
    ax_d.grid(axis='y', color=GRID, lw=0.6)
    ax_d.set_axisbelow(True)
    ax_d.locator_params(axis='y', nbins=3)
    ax_d.tick_params(bottom=False)

    for i, region in enumerate(present):
        pos = df.pos[df.region == region]
        ax_v.vlines(pos, -i - 0.38, -i + 0.38, color=REGIONS[region], lw=0.7)
    ax_v.set_yticks(-np.arange(len(present)),
                    [f'{r}  ({(df.region == r).sum():,})' for r in present], fontsize=11)
    ax_v.set_ylim(-len(present) + 0.5, 0.5)
    ax_v.tick_params(axis='y', length=0)
    ax_v.tick_params(bottom=False)
    ax_v.spines['left'].set_visible(False)
    ax_v.spines['bottom'].set_visible(False)

    draw_gene_models(ax_g, gtf, t.Chromosome, start, end, gene)
    ax_g.set_xlim(start, end)
    digits = 4 if end - start < 30_000 else 3
    ax_g.xaxis.set_major_formatter(lambda x, _: f'{x / 1e6:.{digits}f}')
    ax_g.tick_params(bottom=True)
    ax_g.spines['bottom'].set_visible(True)
    ax_g.set_xlabel(f'{t.Chromosome} position (Mb, hg38)')
    n_coding = (df.region == 'coding').sum()
    fig.suptitle(f'{gene}: {len(df):,} variants ({n_coding:,} coding, {len(df) - n_coding:,} non-coding)',
                 x=0.125, ha='left', y=0.98, style='italic')
    return fig


# --------------------------------------------------------------------------------------
# HTML report
# --------------------------------------------------------------------------------------
CSS = """
body { font-family: system-ui, sans-serif; color: #0b0b0b; background: #fff; max-width: 1100px;
       margin: 2em auto; padding: 0 16px; line-height: 1.45; }
h1 { font-size: 1.6em; } h2 { margin-top: 2em; border-bottom: 1px solid #e1e0d9; }
table { border-collapse: collapse; font-size: 0.9em; margin: 1em 0; }
th, td { padding: 3px 10px; border-bottom: 1px solid #e1e0d9; text-align: right; }
th:first-child, td:first-child { text-align: left; }
img { max-width: 100%; height: auto; margin: 0.5em 0; }
.muted { color: #52514e; }
"""


class Report:
    """Collects HTML sections; figures are saved to ``figdir`` and embedded."""

    def __init__(self, figdir):
        self.figdir = figdir
        self.parts = []

    def h(self, level, text):
        self.parts.append(f'<h{level}>{html.escape(text)}</h{level}>')

    def p(self, text, cls=None):
        attr = f' class="{cls}"' if cls else ''
        self.parts.append(f'<p{attr}>{html.escape(text)}</p>')

    def bullets(self, items):
        self.parts.append('<ul>' + ''.join(f'<li>{html.escape(i)}</li>' for i in items) + '</ul>')

    def table(self, df):
        self.parts.append(df.to_html(border=0, na_rep='', float_format=lambda x: f'{x:,.1f}'))

    def figure(self, fig, name, dpi=150):
        path = self.figdir / f'{name}.png'
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        data = base64.b64encode(path.read_bytes()).decode()
        self.parts.append(f'<img alt="{html.escape(name)}" src="data:image/png;base64,{data}">')

    def write(self, path, title):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            f'<!doctype html>\n<html><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>{html.escape(title)}</title><style>{CSS}</style></head>'
            f'<body>{"".join(self.parts)}</body></html>\n')


def main(argv=None):
    args = parse_args(argv)
    genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    r = Report(args.figdir)

    raw, csq = read_vep_vcf(args.vep_vcf, set(genes))
    variants = raw.drop_duplicates('VariantID').copy()
    variants['type'] = variant_type(variants)
    gtf, mane = load_gtf(args.gtf, genes)
    ann, n_no_gene = annotate(variants, csq, mane)
    add_exon_distance(ann, gtf, mane)
    n = len(ann)
    if n == 0:
        sys.exit(f"No variants on a chromosome with a target gene in {args.vep_vcf}")

    r.h(1, f'VEP report: {args.vep_vcf.name}')
    r.p(f'Input: {args.vep_vcf} (GRCh38). Target genes: {len(genes)}. Gene models: {args.gtf}.', 'muted')
    r.p('Consequences are taken on the MANE Select transcript of a target gene. VEP reports genes within '
        '5 kb; when two target genes qualify, the one the variant lies inside wins, then the nearest. '
        'Variants further than 5 kb from every target gene are "distal" and assigned to the nearest TSS. '
        'Priority when VEP gives several terms: canonical splice site > coding > splice region > UTR > intron.')

    # 1. Variants
    r.h(2, '1. Variants')
    msg = (f'{len(raw):,} VCF records, {len(raw) - len(variants)} duplicate IDs removed '
           f'-> {len(variants):,} unique variants.')
    if n_no_gene:
        msg += f' {n_no_gene:,} variants on chromosomes without a target gene were dropped.'
    r.p(msg)
    r.table(ann.type.value_counts().to_frame('variants'))
    n_other = (ann.region == 'other').sum()
    if n_other:
        r.p(f'{n_other} variants have an unclassified consequence (region "other"); they are not plotted.')
    used = [reg for reg in REGIONS if (ann.region == reg).any()]

    # 2. Genes
    overall, by_type, by_gene = region_tables(ann, used)
    missing = sorted(set(genes) - set(by_gene.index))
    r.h(2, '2. Which genes are covered?')
    r.p(f'{len(by_gene)} of {len(genes)} target genes have variants. Missing: {", ".join(missing) or "none"}.')
    summary = pd.DataFrame({
        'variants': by_gene.sum(axis=1),
        'MANE transcript length (kb)': ((mane.End - mane.Start) / 1e3).round(1),
        'chrom': mane.Chromosome,
    }).loc[by_gene.index]
    summary['variants per kb of gene'] = (summary.variants / summary['MANE transcript length (kb)']).round(1)
    summary['% coding'] = (100 * by_gene.get('coding', 0) / summary.variants).round(1)
    summary['% non-coding'] = (100 - summary['% coding']).round(1)
    r.table(summary)
    r.figure(fig_per_gene(by_gene, used, n), 'variants_per_gene_by_region')

    # 3. Regions
    r.h(2, '3. Where are the variants? Coding vs non-coding')
    r.table(overall)
    n_coding = overall.loc['coding', 'variants'] if 'coding' in overall.index else 0
    r.p(f'coding: {n_coding:,} ({n_coding / n:.1%});  non-coding: {n - n_coding:,} ({1 - n_coding / n:.1%})')
    r.table(by_type)
    detail = detail_counts(ann)
    r.figure(fig_region_and_consequence(overall, detail, 'All variants by region'),
             'variants_by_region_and_consequence')

    r.h(3, 'How deep into introns?')
    r.p('Distance from each intronic or splice-region variant to the nearest MANE exon boundary.')
    d = np.sort(ann.exon_distance.dropna().values)
    if len(d):
        r.figure(fig_exon_distance(d, f'Intronic variants (n={len(d):,}): distance to the nearest exon'),
                 'intronic_distance_to_exon')

    # Summary
    r.h(3, 'Summary')
    top = summary.variants
    items = [
        f'{len(by_gene)} of {len(genes)} target genes have variants; most in '
        + ', '.join(f'{g} ({c:,})' for g, c in top.head(3).items())
        + f'; fewest in {top.index[-1]} ({top.iloc[-1]:,}).',
        f'{n_coding / n:.1%} of variants are coding. Coding consequences: '
        + ', '.join(f'{c:,} {k}' for k, c in detail[CODING_ORDER].items() if c) + '.',
        'Non-coding: ' + ', '.join(f'{c:,} {k}' for k, c in overall.variants.items()
                                   if k != 'coding' and c) + '.',
    ]
    if len(d):
        items.append(f'{(d <= 50).mean():.0%} of intronic / splice-region variants are within 50 bp of an exon.')
    absent = [reg for reg in REGIONS if reg not in used]
    if absent:
        items.append(f'No variants in: {", ".join(absent)}.')
    r.bullets(items)

    # Genome tracks
    r.h(2, '4. Genome tracks per gene')
    r.p('For each gene: variant density, one row of ticks per region class, and the gene models (MANE '
        'transcript of the target gene in black; other protein-coding genes in the window in grey). Tall '
        'boxes are CDS, short boxes UTR exons. Coordinates are hg38.')
    for gene in by_gene.index:
        r.figure(fig_gene_track(ann, gtf, mane, gene), f'tracks/{gene}_variant_track')

    args.annotated.parent.mkdir(parents=True, exist_ok=True)
    ann.to_csv(args.annotated, index=False)
    r.write(args.report_html, f'VEP report {args.vep_vcf.stem}')


if __name__ == '__main__':
    main()
