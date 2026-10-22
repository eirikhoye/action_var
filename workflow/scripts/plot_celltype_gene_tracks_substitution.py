#!/usr/bin/env python3
"""
plot_celltype_gene_tracks_substitution.py
=========================================

Same per-gene, per-cell-type ExPecto tracks as ``plot_celltype_gene_tracks.py``,
but coloured by nucleotide substitution instead of by the sign of the z-score.

Colouring
---------
Substitutions are collapsed to the 6 pyrimidine-reference classes used in
mutational-signature plots (C>A, C>G, C>T, T>A, T>C, T>G; e.g. G>A counts as
C>T). Variants with |z| >= 1 are drawn in full colour, the rest faded, so the
significance threshold stays visible. The legend beside each panel gives the
number of variants per class. Concordant hits are still ringed in black.

Output
------
``<outdir>/<GENE>/<GENE>_<celltype>.<fmt>``

Usage
-----
    python workflow/scripts/plot_celltype_gene_tracks_substitution.py \\
        --genes BRCA1 BRCA2 --formats png svg --threads 16
"""
from plot_celltype_gene_tracks import ROOT, main

if __name__ == "__main__":
    main(description=__doc__,
         outdir=f"{ROOT}/workflow/notebooks/candidate_clin_vars_by_celltype_substitution",
         color_by="substitution")
