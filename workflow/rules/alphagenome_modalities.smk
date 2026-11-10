if 'VEP_ANNOTADED' not in globals():
    raise WorkflowError("alphagenome_modalities needs the annotate_vep module; add it to 'modules'")

AGM = config.get('alphagenome_modalities') or {}
AGM_DIR = AGM.get('outdir', f"{_vep_report_root}_alphagenome_modalities")

