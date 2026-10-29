# rules/sei.smk
# Sei (Troyanskaya lab) sequence-class variant effects for the variants in
# variant_effects, joined back onto that table. Requires the variant_effects module.
#
# One-off setup (see README, "Sei"):
#   git clone https://github.com/FunctionLab/sei-framework && cd sei-framework
#   sh ./download_data.sh            # model + hg19/hg38 FASTA, several GB
#   a Python env with torch, selene-sdk, docopt (plus h5py, pandas)
#
#   sei:
#     dir: /storage/.../sei-framework     # required: the cloned, set-up framework
#     run_prefix: conda run --no-capture-output -n sei   # how to run its Python; "" if active
#     outdir: results/<vep stem>_sei
#     include: noncoding                  # or: all
#     snv_only: false
#     cuda: false                         # true on a GPU node (step 1 is slow on CPU)
#     resources: {}                       # passed to the executor for step 1, e.g.
#                                         #   {slurm_partition: accel, gpus: 1, runtime: 240}

if 'VARIANT_EFFECTS' not in globals():
    raise WorkflowError("sei needs the variant_effects module; add it to 'modules'")

SEI = config.get('sei') or {}
if not SEI.get('dir'):
    raise WorkflowError("sei.dir must point to a set-up sei-framework checkout")
SEI_FRAMEWORK = os.path.abspath(SEI['dir'])
SEI_PREFIX = SEI.get('run_prefix', '')
SEI_OUT = SEI.get('outdir', f"{os.path.splitext(VARIANTS_VEP)[0]}_sei")
SEI_RESOURCES = SEI.get('resources') or {}

SEI_VCF = f"{SEI_OUT}/sei_input.vcf"
_SEI_H5 = f"{SEI_OUT}/chromatin-profiles-hdf5"
SEI_REF = f"{_SEI_H5}/sei_input.ref_predictions.h5"
SEI_ALT = f"{_SEI_H5}/sei_input.alt_predictions.h5"
SEI_LABELS = f"{_SEI_H5}/sei_input_row_labels.txt"
SEI_NPY = f"{SEI_OUT}/sei_input.sequence_class_scores.npy"
SEI_MERGED = f"{SEI_OUT}/variant_effects_sei.tsv"
SEI_LONG = f"{SEI_OUT}/sei_sequence_class_scores_long.tsv"

ALL_OUTPUTS.extend([SEI_MERGED, SEI_LONG])


rule sei_vcf:
    input:
        effects=VARIANT_EFFECTS
    output:
        vcf=SEI_VCF
    params:
        script="workflow/scripts/prepare_sei_vcf.py",
        include=SEI.get('include', 'noncoding'),
        snv_only="--snv-only" if SEI.get('snv_only') else "",
    shell:
        """
        python {params.script} {input.effects} {output.vcf} \
            --include {params.include} {params.snv_only}
        """


# Step 1: 21,907 chromatin-profile predictions for ref and alt (the slow, GPU step).
# The Sei scripts use paths relative to the framework dir, so they are run from there.
rule sei_chromatin:
    input:
        vcf=SEI_VCF
    output:
        ref=SEI_REF,
        alt=SEI_ALT,
        labels=SEI_LABELS,
    params:
        framework=SEI_FRAMEWORK,
        prefix=SEI_PREFIX,
        vcf=lambda wc, input: os.path.abspath(input.vcf),
        outdir=os.path.abspath(SEI_OUT),
        cuda="--cuda" if SEI.get('cuda') else "",
    resources:
        **SEI_RESOURCES
    shell:
        """
        cd {params.framework}
        {params.prefix} sh 1_variant_effect_prediction.sh {params.vcf} hg38 {params.outdir} {params.cuda}
        """


# Step 2: sequence-class scores (alt - ref, histone-normalised); fast, CPU is fine.
# --no-tsv skips Sei's sorted 21,907-column TSVs; merge_sei_scores.py reads the .npy.
rule sei_seqclass:
    input:
        ref=SEI_REF,
        alt=SEI_ALT,
        labels=SEI_LABELS,
    output:
        npy=SEI_NPY
    params:
        framework=SEI_FRAMEWORK,
        prefix=SEI_PREFIX,
        ref=lambda wc, input: os.path.abspath(input.ref),
        alt=lambda wc, input: os.path.abspath(input.alt),
        outdir=os.path.abspath(SEI_OUT),
    shell:
        """
        cd {params.framework}
        {params.prefix} sh 2_varianteffect_sc_score.sh {params.ref} {params.alt} {params.outdir} --no-tsv
        """


rule sei_merge:
    input:
        npy=SEI_NPY,
        labels=SEI_LABELS,
        effects=VARIANT_EFFECTS,
    output:
        merged=SEI_MERGED,
        long=SEI_LONG,
    params:
        script="workflow/scripts/merge_sei_scores.py",
        seqclass=f"{SEI_FRAMEWORK}/model/seqclass.names",
    shell:
        """
        python {params.script} {input.npy} {input.labels} {params.seqclass} \
            {input.effects} {output.merged} --long {output.long}
        """
