# Distal-element ISM with the AlphaGenome Atlas, then Sei

Design for `workflow/rules/alphagenome_ism.smk`. Steps 1–2 run today; steps 3–5 are
stubs to implement.

## Question

Which non-coding positions around each of the 34 BRIDGES genes does AlphaGenome
predict to change that gene's expression in mammary epithelium? Do those positions
form distal elements, and what regulatory activity does Sei assign to them and to
the mutations in them?

The answer gives a map of where variants *could* matter. The observed BRIDGES
variants (from `variant_effects` and `alphagenome`) can then be read against it.

## Pipeline

| Step | Rule | Status | What it does |
|---|---|---|---|
| 1 | `ism_regions` | runs | MANE transcript ± `flank` (default 100 kb) per gene, from the GENCODE feather |
| 2 | `ism_atlas` | runs | Every SNV in the window from the Atlas (`query_interval`), effects on that gene, in 4 kb cached chunks |
| 3 | `ism_elements` | stub | Per position, the max-\|effect\| ALT on the gene's RNA-seq track; keep quantile ≥ 0.99 outside the CDS; merge hits within 50 bp; ≥ 3 positions → element BED + top SNVs |
| 4a | `ism_sei_annotate_elements` | stub | Overlap the elements with Sei's hg38 sequence-class annotation (Zenodo 7113988): what the element *is* (e.g. E6 epithelial enhancer, TF3 FOXA1/AR/ESR1, CTCF). No model run |
| 4b | `ism_sei_variants` | stub | The top ISM SNVs per element through the Sei steps already in `sei.smk`: what those mutations *do* to chromatin |
| 5 | `ism_sei_report` | stub | Element map per gene, AlphaGenome vs Sei concordance, overlap with observed BRIDGES variants |

## Choices and their reasons

- **Atlas, not the model API.** The Atlas returns stored scores, so whole-window ISM is a
  query, not an inference job. That is what `notebooks/alphagenome_atlas_BRCA1.ipynb`
  already does for BRCA1 (120 kb); `ism_atlas` generalises it per gene. The model API
  is only needed for indels.
- **Gene-level RNA-seq effect as the element signal.** It is the quantity that links a
  distal position to *this* gene. AVI and DNase/histone tracks are kept as columns for
  interpretation, not for calling.
- **Calibrated quantile threshold.** The Atlas quantiles make one threshold comparable
  across genes and tracks; raw log-fold changes are not.
- **Sei in two ways.** The annotation overlap (4a) is free and tells you what kind of
  element it is. The variant run (4b) is the like-for-like comparison with AlphaGenome:
  same SNVs, a chromatin-level model trained on different data. Agreement between the two
  is the strongest evidence; disagreement is the interesting part to discuss in Princeton.

## Scale

Three ALTs per position, so 1 Mb of window is about 3 M SNVs. With the default
100 kb flank the 34 windows total roughly 10 Mb (run `ism_regions` to see the exact
number; it prints it). Start with a few genes (`alphagenome_ism.genes`) to time one
4 kb chunk, then scale. The chunk cache makes interrupted runs resumable.

Widening `flank` towards 500 kb catches more distal enhancers but multiplies the query
count. A middle ground is to keep 100 kb for all genes and widen only those whose
step-3 elements cluster near the window edge.

## Open decisions

1. Ontology: mammary epithelial cell (CL:0002327) only, or add breast tissue / MCF-7
   tracks? Check `atlas_client.scorer_metadata()` for which exist per scorer.
2. Element threshold: quantile 0.99 is a starting point. Calibrate it on known elements
   (for example the BRCA1/NBR2 bidirectional promoter) before trusting distal calls.
3. Sei for elements: also run Sei sequence prediction on element BEDs (4 kb centred)
   for elements the genome-wide annotation labels "low signal"?
