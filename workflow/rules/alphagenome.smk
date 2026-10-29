# rules/alphagenome.smk
# Scores the non-coding variants (from variant_effects) with AlphaGenome:
#   - SNVs from the AlphaGenome Atlas (precomputed scores, fast), when backend = atlas
#   - indels, and SNVs the Atlas has no scores for, with the model API
# Requires the variant_effects module. All keys are optional:
#
#   alphagenome:
#     outdir: results/<vep stem>_alphagenome
#     include: noncoding            # or: all
#     exclude_regions: []           # region classes from vep_report, e.g. ["3' UTR"]
#     backend: atlas                # atlas (SNVs) + model (rest); or: model (everything)
#     ontology_terms: [CL:0002327]  # mammary epithelial cell; [] = all biosamples
#     api_key_file: resources/API_KEY   # else $ALPHAGENOME_API_KEY; never commit the key
#     atlas:
#       scorers: [AVI_SCORE, RNA_SEQ, SPLICE_SITES, DNASE, CAGE, CHIP_HISTONE]
#       batch_size: 200
#     model:
#       scorers: [RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, DNASE,
#                 ATAC, CAGE, CHIP_HISTONE, CHIP_TF, POLYADENYLATION]
#       sequence_length: 1MB        # 16KB, 100KB, 500KB or 1MB
#       batch_size: 20
#       all_genes: false            # true: keep gene-level scores for every gene in the window

if 'VARIANT_EFFECTS' not in globals():
    raise WorkflowError("alphagenome needs the variant_effects module; add it to 'modules'")

AG = config.get('alphagenome') or {}
AG_DIR = AG.get('outdir', f"{os.path.splitext(VARIANTS_VEP)[0]}_alphagenome")
AG_BACKEND = AG.get('backend', 'atlas')
if AG_BACKEND not in ('atlas', 'model'):
    raise WorkflowError(f"alphagenome.backend must be 'atlas' or 'model', not {AG_BACKEND!r}")
AG_ONTOLOGY = shlex.quote(','.join(AG.get('ontology_terms', ['CL:0002327'])))
AG_KEY_ARG = f"--api-key-file {shlex.quote(AG['api_key_file'])}" if AG.get('api_key_file') else ""
_AG_ATLAS = AG.get('atlas') or {}
_AG_MODEL = AG.get('model') or {}

AG_SNV = f"{AG_DIR}/input_snv.tsv"
AG_OTHER = f"{AG_DIR}/input_other.tsv"
AG_ATLAS = f"{AG_DIR}/atlas_scores.parquet"
AG_ATLAS_MISSING = f"{AG_DIR}/atlas_missing.tsv"
AG_MODEL = f"{AG_DIR}/model_scores.parquet"

ALL_OUTPUTS.append(AG_MODEL)
if AG_BACKEND == 'atlas':
    ALL_OUTPUTS.append(AG_ATLAS)


rule alphagenome_select:
    input:
        effects=VARIANT_EFFECTS
    output:
        snv=AG_SNV,
        other=AG_OTHER,
    params:
        script="workflow/scripts/select_alphagenome_variants.py",
        include=AG.get('include', 'noncoding'),
        exclude=shlex.quote(','.join(AG.get('exclude_regions') or [])),
    shell:
        """
        python {params.script} {input.effects} {output.snv} {output.other} \
            --include {params.include} --exclude-regions {params.exclude}
        """


rule alphagenome_atlas:
    input:
        snv=AG_SNV
    output:
        scores=AG_ATLAS,
        missing=AG_ATLAS_MISSING,
    params:
        script="workflow/scripts/query_alphagenome_atlas.py",
        scorers=','.join(_AG_ATLAS.get('scorers', ['AVI_SCORE', 'RNA_SEQ', 'SPLICE_SITES',
                                                   'DNASE', 'CAGE', 'CHIP_HISTONE'])),
        batch=_AG_ATLAS.get('batch_size', 200),
        ontology=AG_ONTOLOGY,
        key=AG_KEY_ARG,
    shell:
        """
        python {params.script} variants {input.snv} {output.scores} \
            --missing {output.missing} --scorers {params.scorers} \
            --ontology-terms {params.ontology} --batch-size {params.batch} {params.key}
        """


rule alphagenome_model:
    input:
        # atlas backend: indels + SNVs missing from the Atlas; model backend: everything
        variants=[AG_OTHER, AG_ATLAS_MISSING] if AG_BACKEND == 'atlas' else [AG_SNV, AG_OTHER]
    output:
        scores=AG_MODEL
    params:
        script="workflow/scripts/score_variants_alphagenome.py",
        scorers=','.join(_AG_MODEL.get('scorers', [
            'RNA_SEQ', 'SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS', 'DNASE', 'ATAC',
            'CAGE', 'CHIP_HISTONE', 'CHIP_TF', 'POLYADENYLATION'])),
        length=_AG_MODEL.get('sequence_length', '1MB'),
        batch=_AG_MODEL.get('batch_size', 20),
        all_genes="--all-genes" if _AG_MODEL.get('all_genes') else "",
        ontology=AG_ONTOLOGY,
        key=AG_KEY_ARG,
    shell:
        """
        python {params.script} {output.scores} {input.variants} \
            --scorers {params.scorers} --sequence-length {params.length} \
            --batch-size {params.batch} --ontology-terms {params.ontology} \
            {params.all_genes} {params.key}
        """
