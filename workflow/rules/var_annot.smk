GRCH37_URL = config['GRCh37_url']
GRCH37_path = config['GRCh37_path']
EXPECTO_SCORES = config['expecto_scores']
EXPECTO_RUNFILE = config['expecto_runfile']
ANNOTATED_VARS = config['annotated_variants']
Z_THRESHOLD = config['z_threshold']

GRCH37_path, ANNOTATED_VARS = register_outputs(GRCH37_path, ANNOTATED_VARS)

rule download_gtfs:
    output:
        gtf=GRCH37_path
    params:
        url=GRCH37_URL
    shell:
        """
        wget {params.url} -O {output.gtf}
        """


rule annotate_variants:
    input:
        gtf=GRCH37_path,
        scores=EXPECTO_SCORES,
        run=EXPECTO_RUNFILE
    output:
        csv=ANNOTATED_VARS,
    params:
        script="workflow/scripts/annotate_variants.py",
        z_threshold = Z_THRESHOLD

    shell:
        """
        python {params.script} \
            --runfile {input.run} \
            --zscores {input.scores} \
            --gtf {input.gtf} \
            --out {output.csv} \
            --z-threshold {params.z_threshold}
        """