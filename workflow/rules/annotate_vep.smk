VEP = config.get('vep') or {}

_vep_name = os.path.splitext(os.path.basename(VARIANTS_HG38))[0]
VARIANTS_VEP = VEP.get('output', f"results/{_vep_name}_vep.vcf")

# Input column names, same keys as liftover.columns. Liftover keeps the input's
# column names, so this defaults to liftover.columns when not set.
#
#   vep:
#     columns: {id: SNP, chrom: chrom, pos: snp_chromStart, ref: ref, alt: alt}
VEP_COLUMNS = VEP.get('columns', (config.get('liftover') or {}).get('columns'))
VEP_COLUMN_ARGS = column_args(VEP_COLUMNS, 'vep.columns')

ALL_OUTPUTS.append(VARIANTS_VEP)


rule annotate_vep:
    input:
        var=VARIANTS_HG38
    output:
        vcf=VARIANTS_VEP
    params:
        script="workflow/scripts/annotate_vep.py",
        columns=VEP_COLUMN_ARGS,
    shell:
        """
        python {params.script} {input.var} {output.vcf} {params.columns}
        """


# Report per gene and region class. genes and gtf are required:
#
#   vep:
#     report:
#       output: results/<vep stem>_report.html   # optional
#       genes: [BRCA1, BRCA2, ...]               # target genes (HGNC symbols)
#       gtf: <GENCODE hg38 GTF as feather>       # path or URL
VEP_REPORT = config['vep']['report']
_vep_report_root = os.path.splitext(VEP_REPORT.get('output', os.path.splitext(VARIANTS_VEP)[0] + '_report.html'))[0]
VEP_REPORT_HTML = f"{_vep_report_root}.html"
VEP_ANNOTATED = f"{_vep_report_root}_annotated.csv"
VEP_REPORT_FIGS = f"{_vep_report_root}_figures"

ALL_OUTPUTS.extend([VEP_REPORT_HTML, VEP_ANNOTATED])


rule vep_report:
    input:
        vcf=VARIANTS_VEP
    output:
        html=VEP_REPORT_HTML,
        annotated=VEP_ANNOTATED,
        figures=directory(VEP_REPORT_FIGS)
    params:
        script="workflow/scripts/vep_report.py",
        genes=','.join(VEP_REPORT['genes']),
        gtf=VEP_REPORT['gtf']
    shell:
        """
        python {params.script} {input.vcf} {output.html} \
            --annotated {output.annotated} --figdir {output.figures} \
            --genes {params.genes} --gtf {params.gtf}
        """
