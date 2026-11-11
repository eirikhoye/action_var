#!/usr/bin/env python3
"""
alphagenome_modalities_report.py
================================

HTML report of AlphaGenome effects across modalities for every variant, from the
summary table of ``score_alphagenome_modalities.py`` (one row per variant).

Scale
-----
Every modality is put on a common 0-1 scale, signed for scorers with a direction
(expression, TSS activity, accessibility, histone marks, TF binding). Two choices:

``--scale dataset`` (default)
    Percentile of the variant's |raw score| among all variants in the table, per scorer.
    "Top 1%" then means top 1% of *these* variants, so every modality flags the same
    share overall and the figures show *where* (genes, regions) the strong effects are.
``--scale alphagenome``
    AlphaGenome's own quantile score: rank among common variants genome-wide. For
    gene-level scorers this saturates inside genes (on the 34 BRIDGES genes, 72% of
    SNVs are in the top 1% for splice junctions and 37% for expression), so it does
    not separate variants within a gene panel.

To spread the tail, effects are drawn as ``-log10(1 - |value|)``: 2 = top 1%, 3 = top
0.1%. A variant is "high impact" for a modality when |value| >= ``--high-quantile``.

Modalities
----------
Scorers are grouped into modalities; a variant's modality value is its strongest
scorer in the group. By default only target-biosample tracks are used (breast, see
``--target-ontology`` in the scoring script). Modalities with no target track (3D
contacts, and PRO-cap) use the strongest track in any cell type and are marked "*".

Figures
-------
1. High-impact rate per modality and region class.
2. High-impact variants per gene and modality.
3. Per gene: a multi-modality track. One row per modality, each variant drawn at its
   position with height = effect, coloured by direction; then the region class of each
   variant and the gene models (same layout as the vep_report gene tracks).
4. Per gene: a fingerprint heatmap of the top variants x modalities.

Usage
-----
::

    python alphagenome_modalities_report.py SCORES_PARQUET REPORT_HTML --figdir DIR
        --gtf PATH_OR_URL [--genes G1,G2] [--high-quantile 0.99] [--top-n 25]
"""

import argparse
import html
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vep_report import (AXIS, BLUE, GRID, INK, MUTED, REGIONS, SECONDARY, Report,  # noqa: E402
                        draw_gene_models)

UP, DOWN, UNSIGNED = '#e34948', '#2a78d6', '#eb6834'
MAX_EFFECT = 5.0   # -log10(1 - |q|) is capped here (quantiles are float32)

# (label, scorers, signed)
MODALITIES = [
    ('Expression (RNA-seq)', ['RNA_SEQ'], True),
    ('Splicing', ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS'], False),
    ('Polyadenylation', ['POLYADENYLATION'], False),
    ('TSS activity', ['CAGE', 'PROCAP'], True),
    ('Accessibility', ['DNASE', 'ATAC'], True),
    ('Histone marks', ['CHIP_HISTONE'], True),
    ('TF binding', ['CHIP_TF'], True),
    ('3D contacts', ['CONTACT_MAPS'], False),
]
AVI = ('AVI score', ['AVI_SCORE'], False)
# Modalities whose label lists the scorers actually used (they differ per biosample scope)
SHOW_SCORERS = {'TSS activity', 'Accessibility'}
SHORT = {'CAGE': 'CAGE', 'PROCAP': 'PRO-cap', 'DNASE': 'DNase', 'ATAC': 'ATAC'}
SEQUENTIAL = LinearSegmentedColormap.from_list(
    'blues', ['#ffffff', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
DIVERGING = LinearSegmentedColormap.from_list(
    'div', ['#0d366b', '#2a78d6', '#cde2fb', '#ffffff', '#fbd3d0', '#e34948', '#7d1514'])


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    p.add_argument('scores', type=Path)
    p.add_argument('report_html', type=Path)
    p.add_argument('--figdir', type=Path, required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--genes', default='', help='genes to draw (default: all with variants)')
    p.add_argument('--high-quantile', type=float, default=0.99)
    p.add_argument('--top-n', type=int, default=25)
    p.add_argument('--scope', choices=['target', 'all'], default='target',
                   help="'all': strongest track in any biosample for every modality")
    p.add_argument('--scale', choices=['dataset', 'alphagenome'], default='dataset',
                   help='dataset: percentile among these variants; alphagenome: AlphaGenome quantile')
    return p.parse_args(argv)


# --------------------------------------------------------------------------------------
# Modality table
# --------------------------------------------------------------------------------------
def dataset_percentile(df, qcol):
    """Signed percentile of |raw| among the variants in df (AVI: of the raw score itself)."""
    raw = df[qcol.replace('|q', '|raw')].astype(float)
    if qcol.startswith('AVI_SCORE'):
        return raw.rank(pct=True).to_numpy()
    return (np.sign(raw) * raw.abs().rank(pct=True)).to_numpy()


def modality_table(df, scope, scale='dataset'):
    """Signed score per variant and modality (strongest scorer in the group), on the
    dataset-percentile or AlphaGenome-quantile scale.

    Returns (values DataFrame [variants x modality labels], labels, signed flags)."""
    out, labels, signed = {}, [], {}
    for label, scorers, is_signed in [AVI] + MODALITIES:
        cols, fallback = [], False
        for s in scorers:
            col = f'{s}|q' if scope == 'target' else f'{s}|q_all'
            if col not in df and f'{s}|q' in df:   # scorer without biosamples
                col = f'{s}|q'
            if col in df and df[col].notna().any():
                cols.append(col)
        if not cols:   # no target track at all: any cell type
            cols = [f'{s}|q_all' for s in scorers if f'{s}|q_all' in df]
            fallback = bool(cols)
        if not cols:
            continue
        if scale == 'dataset':
            q = np.column_stack([dataset_percentile(df, c) for c in cols])
        else:
            q = df[cols].to_numpy(dtype=float)
        pick = np.nan_to_num(np.abs(q), nan=-1).argmax(axis=1)
        value = q[np.arange(len(q)), pick]
        name = label
        if label in SHOW_SCORERS:
            name += ' (' + ', '.join(SHORT[c.split('|')[0]] for c in cols) + ')'
        name += ' *' if fallback else ''
        out[name] = value
        labels.append(name)
        signed[name] = is_signed
    return pd.DataFrame(out, index=df.index), labels, signed


def effect(q):
    """-log10(1 - |q|), capped."""
    return np.clip(-np.log10(np.clip(1 - np.abs(q), 10 ** -MAX_EFFECT, 1)), 0, MAX_EFFECT)


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------
def heatmap(ax, values, counts, cmap, vmax, fmt):
    im = ax.imshow(values, aspect='auto', cmap=cmap, vmin=0, vmax=vmax, interpolation='nearest')
    for (i, j), v in np.ndenumerate(values):
        if np.isfinite(v):
            ax.text(j, i, fmt(v, counts[i, j]), ha='center', va='center', fontsize=9,
                    color='white' if v > vmax * 0.6 else INK)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def fig_by_region(mod, df, labels, high):
    regions = [r for r in REGIONS if (df.region == r).any()]
    hit = mod[labels].abs() >= high
    rate = np.array([[hit.loc[df.region == r, m].mean() * 100 for r in regions] for m in labels])
    count = np.array([[hit.loc[df.region == r, m].sum() for r in regions] for m in labels])
    n = [int((df.region == r).sum()) for r in regions]
    fig, ax = plt.subplots(figsize=(1.25 * len(regions) + 3.5, 0.5 * len(labels) + 1.8))
    im = heatmap(ax, rate, count, SEQUENTIAL, max(np.nanmax(rate), 1),
                 lambda v, c: f'{v:.1f}%\n({c:,})' if c else '0')
    ax.set_xticks(range(len(regions)), [f'{r}\n(n={k:,})' for r, k in zip(regions, n)], fontsize=10)
    ax.set_yticks(range(len(labels)), labels, fontsize=11)
    ax.tick_params(length=0)
    ax.xaxis.tick_top()
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(f'% of variants with |score| ≥ {high}', fontsize=10)
    cb.outline.set_visible(False)
    ax.set_title(f'Which modalities do variants in each region affect? (|score| ≥ {high})',
                 pad=40)
    return fig


def fig_by_gene(mod, df, labels, high):
    hit = (mod[labels].abs() >= high)
    counts = hit.groupby(df.gene).sum()
    n = df.gene.value_counts()
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    fig, (ax, ax_n) = plt.subplots(1, 2, figsize=(1.1 * len(labels) + 4, 0.3 * len(counts) + 2.2),
                                   gridspec_kw={'width_ratios': [len(labels), 1.6], 'wspace': 0.04})
    values = counts.values.astype(float)
    im = heatmap(ax, values, counts.values, SEQUENTIAL, max(np.nanquantile(values, 0.98), 1),
                 lambda v, c: f'{c:,}' if c else '')
    ax.set_xticks(range(len(labels)), labels, fontsize=9.5, rotation=35, ha='left')
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(counts)), counts.index, fontsize=10, style='italic')
    ax.tick_params(length=0)
    ax_n.barh(range(len(counts)), n[counts.index], color=AXIS, height=0.7)
    ax_n.set_ylim(len(counts) - 0.5, -0.5)
    ax_n.set_yticks([])
    ax_n.set_xlabel('variants', fontsize=10)
    ax_n.spines['left'].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, location='bottom')
    cb.set_label('number of high-impact variants (grey bars: all variants in the gene)', fontsize=10)
    cb.outline.set_visible(False)
    ax.set_title(f'High-impact variants per gene and modality (|score| ≥ {high})', pad=90)
    return fig


def fig_gene_modalities(df, mod, labels, signed, gtf, mane, gene, high):
    """Multi-modality track: one row per modality, variant effects at their positions."""
    sub, m = df[df.gene == gene], mod[df.gene == gene]
    t = mane.loc[gene]
    lo, hi = min(sub.POS.min(), t.Start), max(sub.POS.max(), t.End)
    pad = int((hi - lo) * 0.03) + 200
    start, end = lo - pad, hi + pad
    thr = effect(high)
    rows = len(labels)
    fig, axes = plt.subplots(
        rows + 2, 1, figsize=(14, 0.62 * rows + 3.4), sharex=True,
        gridspec_kw={'height_ratios': [1] * rows + [0.35, 1.1], 'hspace': 0.12, 'top': 0.84})
    pos = sub.POS.values
    for ax, label in zip(axes, labels):
        q = m[label].values
        e = effect(q)
        strong = np.abs(q) >= high
        ax.scatter(pos[~strong], e[~strong], s=3, color=AXIS, lw=0, rasterized=True)
        if signed[label]:
            color = np.where(q[strong] > 0, UP, DOWN)
        else:
            color = np.full(strong.sum(), UNSIGNED)
        ax.vlines(pos[strong], 0, e[strong], colors=color, lw=0.8, alpha=0.8)
        ax.scatter(pos[strong], e[strong], s=12, c=color, lw=0, zorder=3)
        ax.axhline(thr, color=MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.set_ylim(0, MAX_EFFECT * 1.05)
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color(GRID)
        ax.tick_params(bottom=False)
        ax.set_ylabel(f'{label}  ({int(strong.sum()):,})', rotation=0, ha='right', va='center',
                      fontsize=10.5, color=INK if strong.any() else MUTED)
    # Region class of each variant, one coloured tick per variant
    ax_r = axes[rows]
    ax_r.vlines(pos, 0, 1, colors=sub.region.map(REGIONS).fillna(MUTED).values, lw=0.7)
    ax_r.set_ylim(0, 1)
    ax_r.set_yticks([])
    ax_r.set_ylabel('region', rotation=0, ha='right', va='center', fontsize=10.5)
    for s in ('left', 'bottom'):
        ax_r.spines[s].set_visible(False)
    ax_r.tick_params(bottom=False)
    draw_gene_models(axes[-1], gtf, t.Chromosome, start, end, gene)
    axes[-1].set_xlim(start, end)
    digits = 4 if end - start < 30_000 else 3
    axes[-1].xaxis.set_major_formatter(lambda x, _: f'{x / 1e6:.{digits}f}')
    axes[-1].tick_params(bottom=True)
    axes[-1].spines['bottom'].set_visible(True)
    axes[-1].set_xlabel(f'{t.Chromosome} position (Mb, hg38)')
    # legend: direction + regions present
    present = [r for r in REGIONS if (sub.region == r).any()]
    handles = [plt.Line2D([], [], marker='o', ls='', color=c, ms=6) for c in (UP, DOWN, UNSIGNED)]
    handles += [plt.Rectangle((0, 0), 1, 1, color=REGIONS[r]) for r in present]
    names = ['higher with variant', 'lower with variant', 'change (no direction)'] + present
    fig.legend(handles, names, loc='upper left', bbox_to_anchor=(0.125, 0.925), ncols=6,
               frameon=False, fontsize=9.5, handletextpad=0.3, columnspacing=1.0)
    n_any = int((m[labels].abs() >= high).any(axis=1).sum())
    fig.suptitle(f'{gene}: AlphaGenome effects of {len(sub):,} variants across modalities '
                 f'({n_any:,} high impact in ≥1)', x=0.125, ha='left', y=0.995, va='top', style='italic')
    fig.text(0.125, 0.955, f'Height = -log10(1 - |score|) (dashed: |score| = {high}); '
             'number of high-impact variants per modality in brackets', fontsize=10, color=SECONDARY,
             va='top')
    return fig


def variant_label(r):
    d = r.detail if isinstance(r.detail, str) else (r.region if isinstance(r.region, str) else '')
    return f'{r.CHROM}:{int(r.POS):,} {r.REF[:6]}>{r.ALT[:6]}  {d}'


def fig_fingerprint(df, mod, labels, signed, gene, top_n, high):
    sub, m = df[df.gene == gene], mod.loc[df.gene == gene, labels]
    if sub.empty:
        return None
    order = effect(m.drop(columns=[c for c in m if c.startswith('AVI')], errors='ignore')
                   .to_numpy()).max(axis=1)
    top = np.argsort(-np.nan_to_num(order, nan=-1))[:top_n]
    sub, m = sub.iloc[top], m.iloc[top]
    # colour = direction x effect; unsigned modalities on the positive side
    val = effect(m.to_numpy()) * np.where(m.to_numpy() < 0, -1, 1)
    fig, (ax_r, ax) = plt.subplots(
        1, 2, figsize=(0.75 * len(labels) + 6.5, 0.27 * len(sub) + 2.8),
        gridspec_kw={'width_ratios': [0.25, len(labels)], 'wspace': 0.02,
                     'top': 1 - 1.3 / (0.27 * len(sub) + 2.8)})
    im = ax.imshow(val, aspect='auto', cmap=DIVERGING,
                   norm=TwoSlopeNorm(0, -MAX_EFFECT, MAX_EFFECT), interpolation='nearest')
    hit = np.abs(m.to_numpy()) >= high
    for i, j in zip(*np.nonzero(hit)):
        ax.text(j, i, '•', ha='center', va='center', fontsize=8, color='white')
    ax.set_xticks(range(len(labels)), [l.replace(' (', '\n(') for l in labels], fontsize=8.5,
                  rotation=45, ha='left')
    ax.xaxis.tick_top()
    ax.set_yticks([])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    colors = sub.region.map(REGIONS).fillna(MUTED).tolist()
    ax_r.imshow(np.arange(len(sub))[:, None], aspect='auto',
                cmap=matplotlib.colors.ListedColormap(colors), interpolation='nearest')
    ax_r.set_xticks([])
    ax_r.set_yticks(range(len(sub)), [variant_label(r) for r in sub.itertuples()], fontsize=8.5,
                    family='monospace')
    ax_r.tick_params(length=0)
    for s in ax_r.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label('effect, -log10(1 - |score|), signed by direction\n(splicing, polyA, 3D contacts have no direction: red)', fontsize=9)
    cb.set_ticks([-4, -2, 0, 2, 4], labels=['4 lower', '2', '0', '2', '4 higher'])
    cb.outline.set_visible(False)
    fig.suptitle(f'{gene}: top {len(sub)} variants by strongest modality '
                 f'(• = |score| ≥ {high}; left strip = region class)', fontsize=12, x=0.02, ha='left',
                 y=0.995, va='top', style='italic')
    return fig


# --------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------
def top_table(df, mod, labels, high, n=50):
    mods = [l for l in labels if not l.startswith('AVI')]
    e = effect(mod[mods].to_numpy())
    best = np.nan_to_num(e, nan=-1).argmax(axis=1)
    out = pd.DataFrame({
        'gene': df.gene.values,
        'variant': [f'{r.CHROM}:{int(r.POS)} {r.REF[:10]}>{r.ALT[:10]}' for r in df.itertuples()],
        'region': df.region.values,
        'consequence': df.detail.fillna('').values if 'detail' in df else '',
        'strongest modality': np.array(mods)[best],
        'score': mod[mods].to_numpy()[np.arange(len(df)), best],
        'modalities ≥ threshold': (np.abs(mod[mods]) >= high).sum(axis=1).values,
        'AVI': mod[labels[0]].values if labels[0].startswith('AVI') else np.nan,
    })
    out['_e'] = np.nan_to_num(e, nan=-1).max(axis=1)
    out = out.sort_values(['modalities ≥ threshold', '_e'], ascending=False).drop(columns='_e')
    return out.head(n).reset_index(drop=True)


def main(argv=None):
    args = parse_args(argv)
    from vep_report import load_gtf
    df = pd.read_parquet(args.scores).reset_index(drop=True)
    df['region'] = df.get('region', pd.Series('other', index=df.index)).fillna('other')
    mod, labels, signed = modality_table(df, args.scope, args.scale)
    genes = [g for g in args.genes.split(',') if g] or sorted(df.gene.dropna().unique())
    gtf, mane = load_gtf(args.gtf, genes)
    r = Report(args.figdir)
    high = args.high_quantile
    mods = [l for l in labels if not l.startswith('AVI')]

    r.h(1, 'AlphaGenome effects across modalities')
    r.p(f'Input: {args.scores}. {len(df):,} variants in {df.gene.nunique()} genes '
        f'({(df.backend == "atlas").sum():,} SNVs from the AlphaGenome Atlas, '
        f'{(df.backend == "model").sum():,} indels / MNVs / Atlas misses from the model API).', 'muted')
    if args.scale == 'dataset':
        scale = ('Score = percentile of the variant\'s effect (|raw score|) among all variants in this '
                 'table, per scorer, so "top 1%" means top 1% of these variants and every modality flags '
                 'about the same share overall. AlphaGenome\'s own quantiles (rank among common variants '
                 'genome-wide) saturate inside genes and are kept in the table (columns "|q").')
    else:
        scale = ('Score = AlphaGenome quantile: rank of the variant\'s effect among common variants '
                 'genome-wide. Note that gene-level scores (expression, splicing, polyadenylation) '
                 'saturate inside genes on this scale.')
    r.p(scale + ' Scores are signed for scorers with a direction. Effects are drawn as '
        f'-log10(1 - |score|), so 2 = top 1% and 3 = top 0.1%. "High impact" means |score| ≥ {high}. '
        'Each modality takes the '
        'strongest scorer in its group and, by default, only breast / mammary epithelial tracks; '
        'modalities marked * have no such track and use the strongest track in any cell type. '
        'The AVI score is a composite (AlphaGenome + conservation + protein features) and is shown '
        'separately from the modalities.')
    groups = '; '.join(f'{l}: {", ".join(s)}' for l, s, _ in MODALITIES)
    r.p(f'Modalities (scorers): {groups}.', 'muted')

    hit = mod[mods].abs() >= high
    n_any = int(hit.any(axis=1).sum())
    r.h(2, '1. Overview')
    r.bullets([
        f'{n_any:,} of {len(df):,} variants ({n_any / len(df):.1%}) are high impact in at least one modality; '
        f'{int((hit.sum(axis=1) >= 2).sum()):,} in two or more.',
        'High-impact variants per modality: ' + ', '.join(f'{m} {int(hit[m].sum()):,}' for m in mods) + '.',
    ])
    r.figure(fig_by_region(mod, df, labels, high), 'modality_by_region')
    r.figure(fig_by_gene(mod, df, labels, high), 'modality_by_gene')

    r.h(2, '2. Variants with the broadest effects')
    r.p('Ranked by the number of modalities above the threshold, then by the strongest effect.')
    r.table(top_table(df, mod, labels, high))

    r.h(2, '3. Per gene')
    r.p('Multi-modality track: one row per modality; each variant is drawn at its position with height '
        'equal to its effect. Grey dots are below the threshold; above it, red = higher with the variant, '
        'blue = lower, orange = change without a direction (splicing, polyadenylation, 3D contacts). '
        'Below: region class of each variant and the gene models. The fingerprint heatmap shows the '
        f'top {args.top_n} variants of the gene across all modalities.')
    for gene in genes:
        if not (df.gene == gene).any():
            continue
        r.h(3, gene)
        r.figure(fig_gene_modalities(df, mod, labels, signed, gtf, mane, gene, high),
                 f'tracks/{gene}_modalities')
        fig = fig_fingerprint(df, mod, labels, signed, gene, args.top_n, high)
        if fig is not None:
            r.figure(fig, f'fingerprints/{gene}_top_variants')
        print(f'{gene} done', file=sys.stderr, flush=True)
    r.write(args.report_html, 'AlphaGenome modalities')


if __name__ == '__main__':
    main()