# Download ensembl-vep:
cd resources/
singularity pull ensembl-vep.sif docker://ensemblorg/ensembl-vep:release_111.0

# Test if it works:
singularity exec ensembl-vep.sif vep --help

VEP description:
Ensembl VEP (Variant Effect Predictor) is a tool that takes genomic variants — SNVs, indels, structural variants — and annotates them with their predicted functional consequences.

Core function

Given a variant (e.g. from a VCF file), VEP determines which genes and transcripts it overlaps and predicts what effect it has, such as:

Missense, nonsense, synonymous, frameshift
Splice site disruption (donor/acceptor, splice region)
5'/3' UTR variants
Intronic, intergenic
Regulatory region and non-coding RNA effects
Upstream/downstream gene variants

For coding variants, it gives HGVS notation (c., p.), amino acid changes, and codon changes.

Consequence ranking

VEP assigns a "most severe consequence" per variant using Sequence Ontology terms, and can report all overlapping transcript consequences, not just the worst one — relevant when a variant hits multiple transcript isoforms differently.

Annotation sources it can pull in

This is where it becomes more than a simple overlap tool — VEP supports plugins and custom annotation sources:

Population frequencies (gnomAD, 1000 Genomes)
Clinical significance (ClinVar)
In silico pathogenicity scores (CADD, REVEL, SIFT, PolyPhen)
Conservation scores
Regulatory annotations (from Ensembl Regulatory Build)
Custom BED/BigWig/VCF tracks you supply yourself

Relevant to what you're building with action_var

For non-coding BRCA1/2/PALB2 variants specifically, the parts of VEP that matter most are:

The Regulatory Build annotations (promoter, enhancer, TF binding site overlaps)
Distance-to-transcript metrics for upstream/downstream calls
Custom annotation slots, since you'll likely want to layer in your own ABC-model enhancer-promoter links or Hi-C-derived regulatory elements rather than relying on VEP's default regulatory calls alone
The --pick/--flag_pick type options for collapsing multi-transcript output to something more usable downstream in the pipeline

Since you're running v111 with a RefSeq cache — worth double-checking that RefSeq vs Ensembl transcript sets give consistent regulatory region coordinates, since discrepancies between cache builds are a classic source of the cache-path/output mismatches you mentioned running into.