VEP = config.get('vep') or {}

_vep_name = os.path.splitext(os.path.basename(VARIANTS_HG38))[0]
VARIANTS_VEP = VEP.get('output', f"results/{_vep_name}_vep.vcf")

ALL_OUTPUTS.append(VARIANTS_VEP)


rule annotate_vep:
    input:
        var=VARIANTS_HG38
    output:
        vcf=VARIANTS_VEP
    params:
        script="workflow/scripts/annotate_vep.py"
    shell:
        """
        python {params.script} {input.var} {output.vcf} 
        """