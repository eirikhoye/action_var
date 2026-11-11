if 'VEP_ANNOTATED' not in globals():
    raise WorkflowError("alphagenome_modalities needs the annotate_vep module; add it to 'modules'")

AGM     = config.get('alphagenome_modalities') or {}
AGM_DIR = AGM.get('outdir', f"{_vep_report_root}_alphagenome_modalities")
AGM_TARGET = shlex.quote(','.join(AGM.get('target_ontology', ['CL:0002327', 'UBERON:0008367'])))
AGM_SCORERS = AGM.get('scorers')
AGM_KEY_ARG = f"--api-key-file {shlex.quote(AGM['api_key_file'])}" if AGM.get('api_key_file') else ""
_AGM_ATLAS = AGM.get('atlas') or {}
_AGM_MODEL = AGM.get('model') or {}
_AGM_REPORT = AGM.get('report') or {}

AGM_SNV    = f"{AGM_DIR}/input_snv.tsv"
AGM_OTHER  = f"{AGM_DIR}/input_other.tsv"
AGM_SCRIPT = "workflow/scripts/score_alphagenome_modalities.py"
AGM_ATLAS  = f"{AGM_DIR}/atlas_summary.parquet"
AGM_ATLAS_MISSING = f"{AGM_DIR}/atlas_missing.tsv"
AGM_MODEL = f"{AGM_DIR}/model_summary.parquet"
AGM_SCORES = f"{AGM_DIR}/alphagenome_modality_scores.parquet"
AGM_SCORES_TSV = f"{AGM_DIR}/alphagenome_modality_scores.tsv.gz"
AGM_REPORT_HTML = f"{AGM_DIR}/alphagenome_modalities_report.html"
AGM_REPORT_FIGS = f"{AGM_DIR}/figures"


ALL_OUTPUTS.extend([AGM_SCORES, AGM_REPORT_HTML])


def _agm_scorers(backend):
    if not AGM_SCORERS:
        return ""
    scorers = [s for s in AGM_SCORERS if backend == 'atlas' or s != 'AVI_SCORE']

rule agm_select:
    input:
        annotated=VEP_ANNOTATED
    output:
        snv=AGM_SNV,
        other=AGM_OTHER,
    params:
        regions=shlex.quote(','.join(AGM.get('regions') or [])),
    shell:
        """
        python {AGM_SCRIPT} select {input.annotated} {output.snv} {output.other} \
        --regions {params.regions}
        """


rule agm_atlas:
    input:
        snv=AGM_SNV,
    output:
        scores=AGM_ATLAS,
        missing=AGM_ATLAS_MISSING,
    params:
        scorers=_agm_scorers('atlas'),
        batch=_AGM_ATLAS.get('batch_size', 100),
        workers=_AGM_ATLAS.get('max_workers', 10),
        key=AGM_KEY_ARG,
    shell:
        """
        python {AGM_SCRIPT} score {input.snv} {output.scores} --backend atlas \
            --missing {output.missing} --target-ontology {AGM_TARGET} {params.scorers} \
            --batch-size {params.batch} --max-workers {params.workers} {params.key}
        """


rule agm_model_input:
    # indels / MNVs plus SNVs the Atlas has no scores for, except those whose REF does not
    # match hg38 (the model would score them without checking)
    input:
        other=AGM_OTHER,
        missing=AGM_ATLAS_MISSING,
    output:
        f"{AGM_DIR}/input_model.tsv"
    run:
        import pandas as pd
        missing = pd.read_csv(input.missing, sep='\t')
        if 'atlas_note' in missing:
            bad = missing.atlas_note == 'reference mismatch'
            if bad.any():
                print(f"Skipping {int(bad.sum())} SNVs whose REF does not match hg38 "
                      f"(listed in {input.missing})")
            missing = missing[~bad].drop(columns='atlas_note')
        frames = [pd.read_csv(input.other, sep='\t'), missing]
        pd.concat(frames, ignore_index=True).drop_duplicates('variant_id').to_csv(
            output[0], sep='\t', index=False)

    
rule agm_model:
    input:
        f"{AGM_DIR}/input_model.tsv"
    output:
        scores=AGM_MODEL
    params:
        scorers=_agm_scorers('model'),
        length=_AGM_MODEL.get('sequence_length', '1MB'),
        batch=_AGM_MODEL.get('batch_size', 10),
        workers=_AGM_MODEL.get('max_workers', 5),
        key=AGM_KEY_ARG,
    shell:
        """
        python {AGM_SCRIPT} score {input} {output.scores} --backend model \
            --target-ontology {AGM_TARGET} {params.scorers} --sequence-length {params.length} \
            --batch-size {params.batch} --max-workers {params.workers} {params.key}
        """


rule agm_merge:
    input:
        atlas=AGM_ATLAS,
        model=AGM_MODEL,
        annotated=VEP_ANNOTATED,
    output:
        parquet=AGM_SCORES,
        tsv=AGM_SCORES_TSV,
    shell:
        """
        python {AGM_SCRIPT} merge {output.parquet} {output.tsv} {input.atlas} {input.model} \
        --annotated {input.annotated}
        """


rule agm_report:
    input:
        scores=AGM_SCORES
    output:
        html=AGM_REPORT_HTML,
        figures=directory(AGM_REPORT_FIGS),
    params:
        script="workflow/scripts/alphagenome_modalities_report.py",
        gtf=VEP_REPORT['gtf'],
        genes=shlex.quote(','.join(_AGM_REPORT.get('genes') or [])),
        high=_AGM_REPORT.get('high_quantile', 0.99),
        top_n=_AGM_REPORT.get('top_n', 25),
        scale=_AGM_REPORT.get('scale', 'dataset'),
    shell:
        """
        python {params.script} {input.scores} {output.html} --figdir {output.figures} \
            --gtf {params.gtf} --genes {params.genes} --high-quantile {params.high} \
            --top-n {params.top_n} --scale {params.scale}
        """

