#!/usr/bin/env python3
"""
plot_celltype_gene_tracks.py
============================

Per-cell-type version of the ExPecto gene track plots in
``workflow/notebooks/analyze_variant_annotation.ipynb``.

The notebook plots, for each gene, the mean z-score across the 11 breast
epithelial cell types (``candidate_clin_vars_<GENE>.png``). This script draws
the same track once for every one of the 58 cell types, for every gene.

Input
-----
Same files the notebook uses:

* ``results/expecto_annotated_variants.csv``  (``ann_vars``): one row per
  variant, metadata columns up to ``n_celltypes_abs_z_ge_2``, then one z-score
  column per cell type. ``macro_mean`` is left out (~0 everywhere), leaving 58.
* ``workflow/notebooks/concordant_variants_breast_epithelial.csv``
  (``SELECT_SNPs``): concordant hits, ringed in black.
* ``resources/Homo_sapiens.GRCh37.87.gtf.gz``: gene models.

Output
------
``<outdir>/<GENE>/<GENE>_<celltype>.<fmt>`` for every gene and cell type.

Usage
-----
    python workflow/scripts/plot_celltype_gene_tracks.py \\
        --outdir workflow/notebooks/candidate_clin_vars_by_celltype \\
        --genes BRCA1 BRCA2 --formats png svg --threads 16
"""
import argparse
import os
import re
import subprocess
from functools import lru_cache
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator

ROOT = "/storage/mathelierarea/processed/eirikhoy/action_var"
ANN_VARS = f"{ROOT}/results/expecto_annotated_variants.csv"
CONCORDANT = f"{ROOT}/workflow/notebooks/concordant_variants_breast_epithelial.csv"
GTF = f"{ROOT}/resources/Homo_sapiens.GRCh37.87.gtf.gz"

POS, NEG, NS, INK = "#e34948", "#2a78d6", "#b5b4ae", "#52514e"

# Substitution classes, collapsed to the pyrimidine reference (G>A counts as C>T etc.),
# as in mutational-signature plots. Transitions (C>T, T>C) get the strongest hues.
SUB_COLORS = {"C>A": "#2a78d6", "C>G": "#6250d6", "C>T": "#e34948",
              "T>A": "#eda100", "T>C": "#008300", "T>G": "#e87ba4"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def substitution_class(ref, alt):
    """'C>T'-style class with a pyrimidine reference base."""
    ref, alt = ref.upper(), alt.upper()
    if ref in "GA":
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
    return f"{ref}>{alt}"

# same style as the notebook's figures
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "pdf.fonttype": 42, "ps.fonttype": 42,  # keep text editable in Illustrator
})


@lru_cache(maxsize=None)
def load_gene_models(gene, gtf=GTF):
    """Transcript / exon / CDS rows for one gene from the Ensembl GTF."""
    out = subprocess.run(f"zcat {gtf} | grep -F 'gene_name \"{gene}\";'",
                         shell=True, capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        f = line.split("\t")
        if f[2] not in ("transcript", "exon", "CDS"):
            continue
        attr = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        rows.append({"chrom": f[0], "feature": f[2], "start": int(f[3]), "end": int(f[4]),
                     "strand": f[6], "transcript_id": attr["transcript_id"],
                     "transcript_name": attr.get("transcript_name", attr["transcript_id"]),
                     "biotype": attr.get("transcript_biotype", ""),
                     "ccds": 'tag "CCDS"' in f[8]})
    if not rows:
        raise ValueError(f"{gene} not found in {gtf}")
    return pd.DataFrame(rows)


def plot_gene_track(gene, df, score="mean_z", celltypes=None, highlight=None,
                    window=None, pad=1000, z_thresh=1, n_labels=5,
                    transcripts="protein_coding", max_transcripts=5,
                    symlog=True, gtf=GTF, figsize=(7, 4), color_by="sign"):
    """
    Genome track for one gene: variant z-scores as a lollipop "peak" track
    above the Ensembl gene model. Copied from analyze_variant_annotation.ipynb;
    see the docstring there for the parameters.

    color_by : 'sign' or 'substitution'
        'sign' (as in the notebook): red / blue for positive / negative z,
        grey below z_thresh. 'substitution': colour by substitution class
        (SUB_COLORS), faded below z_thresh, with a legend beside the panel.
    """
    v = df[df["gene_name"] == gene].drop_duplicates("SNP").copy()
    if v.empty:
        raise ValueError(f"No variants for {gene}")
    if celltypes is not None:
        v["z"] = v[list(celltypes)].mean(axis=1)
        ylabel = f"Mean z ({len(celltypes)} cell types)"
    else:
        v["z"] = v[score]
        ylabel = score
    v["pos"] = v["snp_chromEnd"]  # 1-based, matches GTF
    v = v.sort_values("pos")

    lo, hi = window or (v["pos"].min() - pad, v["pos"].max() + pad)
    v = v[v["pos"].between(lo, hi)]
    if v.empty:
        raise ValueError(f"No {gene} variants in window {lo}-{hi}")

    # --- gene model ---
    g = load_gene_models(gene, gtf)
    tx = g[(g["feature"] == "transcript") & (g["end"] >= lo) & (g["start"] <= hi)]
    if isinstance(transcripts, str):
        if transcripts == "protein_coding" and (tx["biotype"] == "protein_coding").any():
            tx = tx[tx["biotype"] == "protein_coding"]
    else:
        tx = tx[tx["transcript_name"].isin(transcripts) | tx["transcript_id"].isin(transcripts)]
    tx = (tx.assign(length=tx["end"] - tx["start"])
            .sort_values(["ccds", "length"], ascending=False)
            .head(max_transcripts))

    n_tx = max(len(tx), 1)
    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw={"height_ratios": [3, 0.35 * n_tx + 0.3], "hspace": 0.05})

    # --- z-score lollipops ---
    sig = v["z"].abs() >= z_thresh
    if color_by == "substitution":
        v["sub"] = [substitution_class(r, a) for r, a in zip(v["ref"], v["alt"])]
        colors = np.array([matplotlib.colors.to_rgba(SUB_COLORS.get(s, NS), 1.0 if strong else 0.25)
                           for s, strong in zip(v["sub"], sig)])
        # draw faded variants first so significant ones sit on top
        order = np.argsort(sig.values, kind="stable")
        v, sig, colors = v.iloc[order], sig.iloc[order], colors[order]
    else:
        colors = np.where(~sig, NS, np.where(v["z"] > 0, POS, NEG))
    ax.vlines(v["pos"], 0, v["z"], colors=colors, lw=0.8, zorder=2)
    ax.scatter(v["pos"], v["z"], c=colors, s=10, lw=0, zorder=3)
    h = v[v["SNP"].isin(set(highlight))] if highlight is not None else v.iloc[:0]
    if len(h):
        ax.scatter(h["pos"], h["z"], s=22, facecolors="none", edgecolors="black",
                   lw=0.6, zorder=4, label=f"concordant hit (n = {len(h)})")
    if color_by == "substitution":
        counts = v["sub"].value_counts()
        handles = [plt.Line2D([], [], marker="o", ls="", color=c, ms=4,
                              label=f"{s}  (n = {counts.get(s, 0):,})") for s, c in SUB_COLORS.items()]
        handles.append(plt.Line2D([], [], marker="o", ls="", color=INK, alpha=0.25, ms=4,
                                  label=f"faded: |z| < {z_thresh}"))
        if len(h):
            handles.append(ax.collections[-1])
        # below the x-axis, where it cannot collide with the variant labels
        fig.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.0),
                   ncol=4, fontsize=7, handletextpad=0.2, columnspacing=1.2)
    elif len(h):
        ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 1), fontsize=7,
                  handletextpad=0.2, borderaxespad=0)
    ax.axhline(0, color=INK, lw=0.5, zorder=1)
    for y in (-z_thresh, z_thresh):
        ax.axhline(y, color=INK, lw=0.5, ls="--", zorder=1)

    ymax = max(v["z"].abs().max(), z_thresh)
    if symlog:
        ax.set_yscale("symlog", linthresh=1)
        dec = 10.0 ** np.arange(0, np.ceil(np.log10(ymax)) + 1)  # 1, 10, 100, ...
        ax.set_yticks(np.r_[-dec[::-1], 0, dec])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}".replace("-", "−")))
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.set_ylim(-ymax * 1.8, ymax * 1.8)
    else:
        ax.set_ylim(-ymax * 1.35, ymax * 1.35)
    ax.set_ylabel(ylabel)
    ax.set_xlim(lo, hi)

    # label the strongest variants: each label goes in the nearest free
    # slot beside its point, alternating up/down, kept inside the panel
    fig.canvas.draw()
    to_px, to_data = ax.transData.transform, ax.transData.inverted().transform
    line_px = 9 * fig.dpi / 72  # one label height
    near_px = 90 * fig.dpi / 72  # labels closer than this in x may collide
    y0, y1 = ax.bbox.y0 + line_px / 2, ax.bbox.y1 - line_px / 2
    placed = []
    for _, r in v.loc[v["z"].abs().nlargest(n_labels).index].iterrows():
        px, py = to_px((r["pos"], r["z"]))
        for k in (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6):
            ty = min(max(py + k * line_px, y0), y1)
            if not any(abs(qx - px) < near_px and abs(qy - ty) < line_px for qx, qy in placed):
                break
        placed.append((px, ty))
        ax.annotate(r["SNP"], (r["pos"], r["z"]), xytext=to_data((px + 12, ty)),
                    fontsize=6, color=INK, va="center",
                    arrowprops=dict(arrowstyle="-", color=INK, lw=0.4, shrinkA=0, shrinkB=2))

    ax.set_title(f"{gene}  ({len(v):,} variants, {int(sig.sum())} with |z| ≥ {z_thresh})",
                 loc="left", fontsize=9)
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False)

    # TSS used by ExPecto, drawn through both panels
    for tss, cage in v[["tss", "cage_used"]].drop_duplicates().itertuples(index=False):
        for a in (ax, axg):
            a.axvline(tss, color=INK, lw=0.6, ls=":", zorder=0)
        axg.annotate("TSS (CAGE)" if cage else "TSS", (tss, 1), xycoords=("data", "axes fraction"),
                     xytext=(2, 0), textcoords="offset points", fontsize=6, color=INK, va="top")

    # --- transcripts ---
    for i, t in enumerate(tx.itertuples()):
        y = -i
        parts = g[g["transcript_id"] == t.transcript_id]
        axg.hlines(y, max(t.start, lo), min(t.end, hi), color=INK, lw=0.6)
        # strand chevrons along the intron line
        step = (hi - lo) / 40
        for x in np.arange(max(t.start, lo) + step / 2, min(t.end, hi), step):
            axg.text(x, y, ">" if t.strand == "+" else "<", fontsize=5, color=INK,
                     ha="center", va="center")
        for p in parts.itertuples():
            if p.feature == "transcript":
                continue
            hgt = 0.35 if p.feature == "CDS" else 0.18  # UTR/non-coding exons thinner
            axg.fill_between([p.start, p.end], y - hgt, y + hgt, color=INK, lw=0, zorder=3)
        axg.text(lo, y + 0.3, t.transcript_name, fontsize=6, color=INK, va="bottom")
    axg.set_ylim(-n_tx + 0.4, 1.0)
    axg.set_yticks([])
    axg.spines[["top", "right", "left"]].set_visible(False)

    axg.set_xlim(lo, hi)
    axg.xaxis.set_major_locator(MaxNLocator(5))
    axg.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    axg.set_xlabel(f"chr{v['chrom'].iloc[0]} position (GRCh37)")
    return fig, (ax, axg)


def load_inputs(ann_vars_path=ANN_VARS, concordant_path=CONCORDANT):
    """ann_vars, the 58 cell-type columns and SELECT_SNPs, as in the notebook."""
    ann_vars = pd.read_csv(ann_vars_path)
    meta_end = ann_vars.columns.get_loc("n_celltypes_abs_z_ge_2")
    cell_cols = [c for c in ann_vars.columns[meta_end + 1:] if c != "macro_mean"]  # macro_mean is ~0 everywhere
    select_snps = pd.read_csv(concordant_path)["SNP"]
    return ann_vars, cell_cols, select_snps


def plot_gene(job):
    """All cell-type tracks for one gene (one worker process per gene)."""
    gene, ann_vars, cell_cols, select_snps, outdir, formats, dpi, color_by = job
    gene_dir = os.path.join(outdir, gene)
    os.makedirs(gene_dir, exist_ok=True)
    for cell in cell_cols:
        fig, (ax, _) = plot_gene_track(gene, ann_vars, score=cell, highlight=select_snps, color_by=color_by)
        ax.set_title(f"{cell}:  {ax.get_title(loc='left')}", loc="left", fontsize=9)
        ax.set_ylabel(f"z ({cell})")
        for fmt in formats:
            fig.savefig(os.path.join(gene_dir, f"{gene}_{cell}.{fmt}"), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
    return gene, len(cell_cols)


def main(description=__doc__, outdir=f"{ROOT}/workflow/notebooks/candidate_clin_vars_by_celltype",
         color_by="sign"):
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=outdir)
    parser.add_argument("--ann-vars", default=ANN_VARS)
    parser.add_argument("--concordant", default=CONCORDANT)
    parser.add_argument("--genes", nargs="+", help="genes to plot (default: all genes in ann_vars)")
    parser.add_argument("--formats", nargs="+", default=["png"], help="e.g. png svg pdf")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    ann_vars, cell_cols, select_snps = load_inputs(args.ann_vars, args.concordant)
    genes = args.genes or list(ann_vars["gene_name"].unique())
    print(f"{len(genes)} genes x {len(cell_cols)} cell types -> {args.outdir}", flush=True)

    jobs = [(g, ann_vars[ann_vars["gene_name"] == g], cell_cols, select_snps,
             args.outdir, args.formats, args.dpi, color_by) for g in genes]
    with Pool(args.threads) as pool:
        for gene, n in pool.imap_unordered(plot_gene, jobs):
            print(f"  {gene}: {n} plots", flush=True)


if __name__ == "__main__":
    main()
