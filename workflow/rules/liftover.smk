# rules/liftover.smk
# Lifts VARIANTS between genome builds (default hg19 -> hg38).
# All globals below are defined here, so they only exist when "liftover" is
# listed in config['modules']. Every config key is optional:
#
#   liftover:
#     from_build: hg19
#     to_build: hg38
#     chain: resources/hg19ToHg38.over.chain.gz   # omit to let pyliftover download it
#     output: results/variants_hg38.csv           # omit to derive from VARIANTS

LIFTOVER = config.get('liftover', {})
LIFTOVER_FROM = LIFTOVER.get('from_build', VARIANTS_BUILD)
LIFTOVER_TO = 'hg38'
LIFTOVER_CHAIN = LIFTOVER.get('chain')

if LIFTOVER_FROM != LIFTOVER_TO:
    _variants_root, _variants_ext = os.path.splitext(VARIANTS)
    VARIANTS_HG38 = LIFTOVER.get('output', f"{_variants_root}_{LIFTOVER_TO}{_variants_ext}")
    ALL_OUTPUTS.append(VARIANTS_HG38)

    rule liftover:
        input:
            variants=VARIANTS
        output:
            lifted=VARIANTS_HG38
        params:
            script="workflow/scripts/genome_liftover.py",
            from_build=LIFTOVER_FROM,
            to_build=LIFTOVER_TO,
            chain=f"--chain {LIFTOVER_CHAIN}" if LIFTOVER_CHAIN else ""
        shell:
            """
            python {params.script} {input.variants} {output.lifted} \
                --from-build {params.from_build} --to-build {params.to_build} {params.chain}
            """