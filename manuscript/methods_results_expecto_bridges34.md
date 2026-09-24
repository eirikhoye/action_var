# Methods

## Variant set and cell-type-resolved effect predictions

We analysed [N = 10,066] single-nucleotide variants (SNVs) in the regulatory regions of the 34 breast cancer susceptibility genes from the BRIDGES panel [ref: Dorling et al., 2021]. [Describe how the variants were chosen, e.g. source database / saturation / population frequency filter.] All coordinates are on GRCh37/hg19. We assigned each variant to the gene whose transcription start site (TSS) was closest. TSSs were taken from CAGE data when available ([8,663] variants) and from the gene annotation otherwise ([1,403] variants). All variants lay within 20 kb of their assigned TSS (median distance 11.3 kb). This nearest-TSS rule assigned variants to 55 genes: the BRIDGES genes plus 27 neighbouring genes whose TSS was closer to some variants than the BRIDGES gene's TSS. After removing duplicate rows, one variant (1_45805963_C_T) was equidistant (127 bp) to the TSSs of *TOE1* and *MUTYH*; we kept the first entry after sorting by distance (*TOE1*).

For each variant we obtained predicted effects on expression in 58 cell types from a normal human breast single-cell atlas [ref]. The predictions came from [model: ExPecto-sc / Seqweaver-based; ref] and are given as z-scores, [describe the scaling, e.g. standardised against a null distribution of variant effects]. The 58 cell types comprise 11 epithelial (luminal secretory, luminal hormone-responsive and basal), 13 stromal/vascular (fibroblast, vascular and lymphatic endothelial, pericyte and vascular smooth muscle), 14 myeloid and 20 lymphoid populations. We excluded an aggregated macrophage column (`macro_mean`) because it summarises other columns.

## Principal component analysis of cell types

We stored the z-scores as a cell type × variant matrix (58 × 10,066) with AnnData (v[x]). We ran principal component analysis (PCA) with Scanpy (v[x]; `sc.tl.pca`, ARPACK solver, 10 components). Each variant was mean-centred across cell types but not scaled to unit variance, so variants with large predicted effects contribute more. We show PC1 and PC2 as a biplot. Cell types are plotted by their PC scores. The 15 variants with the largest loading magnitude (Euclidean norm of the PC1 and PC2 loadings) are drawn as arrows, rescaled to the range of the cell-type scores.

We took the 30 variants with the largest PC1/PC2 loading magnitude and plotted their z-scores across all 58 cell types as a heatmap. Rows (variants) and columns (cell types) were ordered by hierarchical clustering (seaborn `clustermap`; Euclidean distance, average linkage).

## Gene-level summary of variants with large predicted effects

We converted the matrix to long format, one row per variant–cell type pair (583,828 pairs), and kept pairs with |z| > 1. Genes were labelled by their expected role in breast carcinogenesis. The 27 BRIDGES genes that are tumour suppressors (loss-of-function risk genes) were labelled "suppressor". *AKT1* was labelled "oncogene". All other genes, including the neighbouring genes picked up by the nearest-TSS rule, were labelled "unclassified". For suppressor genes we defined the expected direction as a decrease in expression (z < 0). For *AKT1* we defined it as an increase (z > 0). A variant–cell type pair was called *directionally concordant* when the sign of its z-score matched the expected direction.

For each gene we report the following, over pairs with |z| > 1:

- the number of distinct variants;
- the genomic span of all tested variants assigned to the gene;
- the lead variant, i.e. the one with the largest |z| in any cell type;
- the largest number of cell types in which any single variant has |z| > 1;
- the mean percentage of directionally concordant cell types per variant.

Last, we selected variants that were directionally concordant in at least one of the 11 breast epithelial cell types, which are the likely cells of origin of breast carcinoma. Unclassified genes were excluded from this step. For each selected variant we recorded how many epithelial cell types met the criterion, its largest |z| and its mean z in those cell types. Analyses used Python 3.12 with pandas, NumPy, AnnData, Scanpy, matplotlib and seaborn.

# Results

## Predicted variant effects are sparse and heavy-tailed

Across the 583,828 variant–cell type pairs, predicted effects were centred on zero (mean z = −0.002, median 0.005, s.d. 1.43) but had heavy tails (range −36.6 to 68.1). Only 8.95% of pairs had |z| > 1, 4.79% had |z| > 2 and 0.30% had |z| > 10. Pairs with |z| > 1 came from 1,895 of the 10,066 variants (18.8%), and pairs with |z| > 2 from 1,048 variants (10.4%). Most variants were therefore predicted to have little effect in any cell type, and a small subset had large effects.

## Cell types separate by lineage along the main axes of variant-effect variation

PC1 and PC2 explained 44.4% and 20.8% of the variance in predicted variant effects, 65.2% combined (Figure [X]A). PC1 separated lymphoid cells from epithelial and stromal cells. T-cell, NK-cell and B-cell populations had the highest PC1 scores (16 to 57; T_prol 57.1, CD8_activated 57.0). Fibroblasts, endothelial cells and epithelial cells had the lowest (−14 to −40). PC2 separated myeloid cells from luminal epithelial cells. Neutrophils (51.5), monocytes and macrophage subsets (34 to 41) scored highest on PC2. Luminal secretory and hormone-responsive epithelial cells scored lowest (−47 to −21; Lumsec_prol −46.8), and stromal cells sat near zero. The main differences in predicted variant effects between cell types therefore follow the major cell lineages of the breast: epithelial, stromal, myeloid and lymphoid.

The variants that contributed most to PC1 and PC2 came from a few genes. Among the 30 top-loading variants, 10 were assigned to *FANCM*, 6 to *PTEN*, 5 to *MSH6*, 2 each to *STK11* and *ATM*, and 1 each to *TOE1*, *EPM2AIP1*, *NDUFAB1*, *GEN1* and *CDH1* (Figure [X]B). The strongest was 10_89624223_A_G in *PTEN* (mean z = 45.6; max z = 68.1). Like the other top *PTEN* variants, it was predicted to increase expression in every cell type, most strongly in lymphoid cells (mean z 57.5 in lymphoid vs. 38.4 in epithelial). The top *FANCM* variants had large effects in both directions (e.g. 14_45605454_G_A, mean z = −15.5; 14_45605532_C_G, mean z = 14.9). *STK11* variant 19_1207160_A_G and *MSH6* variant 2_48026020_C_T loaded mainly on PC2 and had their largest effects in myeloid cells (mean z 13.1 and −10.8, respectively).

## Variants with large predicted effects cluster in DNA-repair and tumour-suppressor genes

Of the 1,895 variants with |z| > 1 in at least one cell type, 1,512 (79.8%) were assigned to 26 suppressor genes, 8 to *AKT1* and 375 to 22 unclassified genes (Table [Y]; Figure [X]C). Suppressor genes made up 63.0% of all tested variants (6,340 of 10,066), so variants with large predicted effects were over-represented among them. The genes with the most such variants were *MSH2* (220 variants; lead variant 2_47630103_C_T, |z| = 26.1), *MSH6* (159), *PALB2* (146), *FANCM* (115), *MLH1* (106) and *PMS2* (97), which include all four mismatch-repair genes. The largest single effects were in *PTEN* (10_89624223_A_G, |z| = 68.1), *FANCM* (14_45605454_G_A, |z| = 36.6), *EPM2AIP1* (3_37034760_A_G, |z| = 28.2), *MSH2* (|z| = 26.1), *ATM* (11_108093617_T_C, |z| = 25.0) and *MLH1* (3_37034997_C_T, |z| = 24.0). *EPM2AIP1* shares a bidirectional promoter with *MLH1*. For 25 of the 26 suppressor genes with any variant above threshold, at least one variant had |z| > 1 in 52–58 cell types, so strong predicted effects were usually not confined to one cell type.

## Direction of predicted effects is not enriched for the expected loss-of-function direction

Across suppressor genes, 53.2% of variant–cell type pairs with |z| > 1 predicted reduced expression. This is close to the 52.5% negative fraction among all pairs with |z| > 1, whatever the gene's role. Per-gene mean concordance ranged from 22% (*TP53*) to 74% (*GEN1*), with a median of 49%. Predicted effects were therefore not biased towards the loss-of-expression direction expected for risk variants in tumour-suppressor genes. Instead, each gene had a mix of variants predicted to increase or decrease expression.

## Candidate expression-reducing variants in breast epithelial cells

Restricting to the 11 breast epithelial cell types, 633 variants in 25 genes were directionally concordant in at least one epithelial cell type (5,085 variant–cell type pairs). These effects were usually shared across the epithelial compartment. The median variant was concordant in all 11 epithelial cell types (329 of 633 variants), and 453 variants (71.6%) were concordant in at least 6. *MSH2* (112 variants), *PALB2* (70), *PMS2* (62), *FANCM* (57) and *MLH1* (50) contributed the most candidates. *BRCA1* contributed 27 and *BRCA2* 10. The strongest candidates were 3_37034997_C_T in *MLH1* (max |z| = 23.1 in epithelial cells), 10_89624335_C_T in *PTEN* (21.6), 14_45605454_G_A in *FANCM* (20.7), 2_47630215_G_T in *MSH2* (16.1) and 17_41277374_G_A in *BRCA1* (14.4). Each was predicted to reduce expression in all 11 epithelial cell types. These variants are candidate regulatory risk alleles for functional follow-up.

# Figure and table legends (draft)

**Figure [X]. Cell-type-resolved predicted effects of regulatory variants in BRIDGES breast cancer genes.** (A) PCA biplot of 58 breast cell types (points) based on predicted expression effects (z-scores) of 10,066 variants. Arrows show the 15 variants with the largest PC1/PC2 loadings, labelled with their assigned gene. (B) Hierarchically clustered heatmap of z-scores for the 30 variants with the largest PC1/PC2 loading magnitude across all 58 cell types. (C) Number of variants with |z| > 1 in at least one cell type, per gene. Bar colour shows the percentage of directionally concordant cell types (see Methods). Annotations give the lead variant's |z| and the largest number of cell types in which any single variant has |z| > 1.

**Table [Y]. Gene-level summary of variants with |z| > 1** (from `report_sig_variants_per_gene.csv`). Columns: gene role, number of variants, genomic span, lead variant, lead variant max |z|, largest number of cell types per variant, percentage concordant.
