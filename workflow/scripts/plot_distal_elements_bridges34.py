#!/usr/bin/env python3
"""
plot_distal_elements_bridges34.py
=================================

Figure 3b of ``workflow/notebooks/alphagenome_seminar_figures.ipynb``
("Beyond the gene: distal regulatory elements") for every BRIDGES 34 gene.

The plotting code is not copied: the setup cells, helper functions and
``figure_distal_ism`` are read from the notebook and executed here, so the
figures always match the notebook. Only imports, constants and function
definitions are taken from the helper cells; cells that make the notebook's
own figures are skipped.

For each gene the model predicts DNase/H3K27ac over the 1 Mb window around the
MANE TSS, the 60 strongest DNase peaks are called, and each is mutagenised in
the AlphaGenome Atlas (400 bp, cached under ``data/alphagenome_atlas/<GENE>/``).
Atlas queries for several genes are run in parallel first (``--threads``); the
figures are then drawn one by one.

Output
------
``<outdir>/fig03b_atlas_distal_elements_<GENE>.png`` / ``.pdf``
``<outdir>/distal_elements_<GENE>.csv``  (one row per element, as ``DISTAL``)

Usage
-----
    python workflow/scripts/plot_distal_elements_bridges34.py
    python workflow/scripts/plot_distal_elements_bridges34.py --genes BRCA2 TP53 --threads 2
"""
import argparse
import ast
import concurrent.futures
import json
import os
import pathlib
import traceback

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = pathlib.Path("/storage/mathelierarea/processed/eirikhoy/action_var")
NOTEBOOK = ROOT / "workflow/notebooks/alphagenome_seminar_figures.ipynb"
OUTDIR = ROOT / "workflow/notebooks/seminar_figures/distal_elements_bridges34"

# Cells run as a whole (setup, style, gene annotation, track drawing helpers),
# found by a line they contain. Bare expressions are dropped, except rcParams.
FULL_CELLS = ["API_KEY =", "SLIDE =", "BRIDGES = MANE", "def style_track"]
# Cells from which only imports, function definitions and these constants are taken.
DEF_CELLS = ["def get_track", "def fetch_region", "def figure_distal_ism"]
CONSTANTS = {"MODEL_OUTPUTS", "MODEL_TRACKS", "ATLAS_SCORERS", "CHUNK", "AVI", "RNA",
             "SPLICE", "SIGNED", "REGION_ORDER", "DNASE"}


def notebook_code(path=NOTEBOOK):
    """Source of the notebook statements needed for figure_distal_ism."""
    cells = ["".join(c["source"]) for c in json.load(open(path))["cells"] if c["cell_type"] == "code"]

    def find(marker):
        hits = [s for s in cells if marker in s]
        if not hits:
            raise ValueError(f"no notebook cell contains {marker!r}")
        return hits[0]

    keep = []
    for marker in FULL_CELLS:
        src = find(marker)
        for node in ast.parse(src).body:
            seg = ast.get_source_segment(src, node)
            if isinstance(node, ast.Expr) and not seg.startswith("plt.rcParams"):
                continue
            keep.append(seg)
    for marker in DEF_CELLS:
        src = find(marker)
        for node in ast.parse(src).body:
            names = {t.id for n in ([node] if isinstance(node, ast.Assign) else [])
                     for t in ast.walk(n.targets[0]) if isinstance(t, ast.Name)}
            if (isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))
                    or (names and names <= CONSTANTS)):
                keep.append(ast.get_source_segment(src, node))
    return "\n\n".join(keep)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=str(OUTDIR))
    parser.add_argument("--genes", nargs="+", help="genes to plot (default: all 34 BRIDGES genes)")
    parser.add_argument("--threads", type=int, default=4, help="genes queried in parallel")
    parser.add_argument("--n-elements", type=int, default=60)
    args = parser.parse_args()

    # The notebook uses paths relative to workflow/notebooks.
    os.chdir(NOTEBOOK.parent)
    ns = {"display": print}
    exec(compile(notebook_code(), str(NOTEBOOK), "exec"), ns)
    plt, genome, dna_client = ns["plt"], ns["genome"], ns["dna_client"]

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    genes = args.genes or ns["BRIDGES_GENES"]

    def prefetch(gene):
        """Warm the Atlas cache for one gene (same calls as figure_distal_ism)."""
        t = ns["BRIDGES"].loc[gene]
        window = genome.Interval(t.Chromosome, t.tss, t.tss).resize(dna_client.SEQUENCE_LENGTH_1MB)
        output = ns["dna_model"].predict_interval(
            interval=window, ontology_terms=[ns["ONTOLOGY"]],
            requested_outputs={dna_client.OutputType.DNASE, dna_client.OutputType.CHIP_HISTONE},
        )
        ns["element_ism"](ns["call_elements"](output, window, args.n_elements), gene, t.Chromosome)
        return gene

    print(f"Atlas queries for {len(genes)} genes ({args.threads} in parallel)", flush=True)
    failed = {}
    with concurrent.futures.ThreadPoolExecutor(args.threads) as pool:
        futures = {pool.submit(prefetch, g): g for g in genes}
        for f in concurrent.futures.as_completed(futures):
            gene = futures[f]
            try:
                f.result()
                print(f"  cached {gene}", flush=True)
            except Exception:
                failed[gene] = traceback.format_exc()
                print(f"  FAILED {gene}\n{failed[gene]}", flush=True)

    for gene in genes:
        if gene in failed:
            continue
        try:
            fig, elements = ns["figure_distal_ism"](gene, n_elements=args.n_elements)
        except Exception:
            failed[gene] = traceback.format_exc()
            print(f"  FAILED {gene}\n{failed[gene]}", flush=True)
            continue
        for ext in ("png", "pdf"):
            fig.savefig(outdir / f"fig03b_atlas_distal_elements_{gene}.{ext}", dpi=250, bbox_inches="tight")
        plt.close(fig)
        elements.sort_values("summit").to_csv(outdir / f"distal_elements_{gene}.csv", index=False)
        print(f"  plotted {gene}", flush=True)

    print(f"done: {len(genes) - len(failed)} of {len(genes)} genes -> {outdir}")
    if failed:
        print("failed:", ", ".join(failed))


if __name__ == "__main__":
    main()
