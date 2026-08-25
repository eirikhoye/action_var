# /rules/alphagenome.smk

def alphagenome_outputs():
    return [SCORES]

ALL_OUTPUTS.extend(alphagenome_outputs())



rule score_variants_alphagenome:
    input:
        tsv=VARIANTS
    output:
        scores=SCORES
    params:
        script="workflow/scripts/score_variants_alphagenome.py"
    shell:
        """
        python {params.script} {input.tsv} {output.scores}
        """


