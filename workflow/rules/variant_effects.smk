# rules/variant_effects.smk
# Applies the variant-effect criteria from Denise's paper (LoF, damaging missense,
# splicing, optional ExPectoSC expression) to the annotate_vep output.
# Requires the annotate_vep module. All keys are optional:
#
#   variant_effects:
#     output: results/<vep stem>_effects.tsv
#     summary: results/<vep stem>_effects_summary.tsv
#     thresholds:                     # defaults = Denise's table
#       polyphen_min: 0.85
#       revel_min: 0.6
#       metarnn_min: 0.8
#       dbscsnv_min: 0.9              # both ada and rf
#       maxentscan_abs_min: 2
#       spliceai_min: 0.5
#       expectosc_abs_z: 1
#       missense_min_tools: 1         # tools that must agree for damaging_missense
#       splice_min_tools: 1
#     expectosc:                      # optional: z-scores from an ExPectoSC run
#       table: .../ExpectoSC_bridges__breast.csv
#       id_col: SNP                   # must match the VEP VCF IDs (= input VariantID)
#       celltypes: []                 # empty = all numeric columns
import json

if 'VARIANTS_VEP' not in globals():
    raise WorkflowError("variant_effects needs the annotate_vep module; add it to 'modules'")

EFFECTS = config.get('variant_effects') or {}
_vep_root = os.path.splitext(VARIANTS_VEP)[0]
VARIANT_EFFECTS = EFFECTS.get('output', f"{_vep_root}_effects.tsv")
VARIANT_EFFECTS_SUMMARY = EFFECTS.get('summary', f"{_vep_root}_effects_summary.tsv")
_EXPECTOSC = EFFECTS.get('expectosc') or {}

ALL_OUTPUTS.extend([VARIANT_EFFECTS, VARIANT_EFFECTS_SUMMARY])


rule classify_variant_effects:
    input:
        vcf=VARIANTS_VEP,
        regions=VEP_ANNOTATED,
        expectosc=[_EXPECTOSC['table']] if _EXPECTOSC.get('table') else [],
    output:
        tsv=VARIANT_EFFECTS,
        summary=VARIANT_EFFECTS_SUMMARY,
    params:
        script="workflow/scripts/classify_variant_effects.py",
        genes=','.join(VEP_REPORT['genes']),
        thresholds=shlex.quote(json.dumps(EFFECTS.get('thresholds') or {})),
        expectosc=(f"--expectosc {shlex.quote(_EXPECTOSC['table'])} "
                   f"--expectosc-id-col {shlex.quote(_EXPECTOSC.get('id_col', 'SNP'))} "
                   + (f"--expectosc-celltypes {shlex.quote(','.join(_EXPECTOSC['celltypes']))}"
                      if _EXPECTOSC.get('celltypes') else "")) if _EXPECTOSC.get('table') else "",
    shell:
        """
        python {params.script} {input.vcf} {output.tsv} --genes {params.genes} \
            --regions {input.regions} --thresholds {params.thresholds} \
            {params.expectosc} --summary {output.summary}
        """
