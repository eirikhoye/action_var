INPUT  = "/storage/mathelierarea/processed/eirikhoy/action_var/data/BRCA1_BRCA2_PALB2_varianter_2025_02_10.csv"
VCF    = "/storage/mathelierarea/processed/eirikhoy/action_var/results/clinical_gene_variants.vcf"
TSV = VCF.removesuffix('.vcf') + "_alphagenome.tsv"
SCORES = "/storage/mathelierarea/processed/eirikhoy/action_var/results/clinical_gene_scores.csv"

def alphagenome_outputs():
    return TSV, SCORES


rule convert_hgvs2vcf:
    input:
        txt=INPUT
    
    output:
        vcf=VCF,
        tsv=TSV

    params:
        prefix=lambda wc, output: output.vcf.removesuffix(".vcf"),
        script="workflow/scripts/hgvs2vcf.py"

    shell:
        """
        python {params.script} {input.txt} -o {params.prefix}
        """


rule score_variants_alphagenome:
    input:
        tsv=TSV
    output:
        scores=SCORES
    params:
        script="workflow/scripts/score_variants_alphagenome.py"
    shell:
        """
        python {params.script} {input.tsv} {output.scores}
        """


