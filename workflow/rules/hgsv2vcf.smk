import os
base, ext = os.path.splitext(VARIANTS)

PREFIX = base + "_converted"
VCF = PREFIX + ".vcf"
TSV = PREFIX + "_alphagenome.tsv"

VCF, TSV = register_outputs(VCF, TSV)

rule convert_hgvs2vcf:
    input:
        txt=VARIANTS
    
    output:
        vcf=VCF,
        tsv=TSV

    params:
        prefix=PREFIX,
        script="workflow/scripts/hgvs2vcf.py"

    shell:
        """
        python {params.script} {input.txt} -o {params.prefix}
        """

VARIANTS = TSV