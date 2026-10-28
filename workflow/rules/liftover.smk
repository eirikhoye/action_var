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
#     columns:                                    # omit any key to use the default name
#       id: VariantID
#       chrom: Chromosome
#       pos: Position
#       ref: Ref
#       alt: Alt
#     extra_pos_columns:                          # other coordinate columns to lift,
#       snp_chromEnd:                             # as <position col>: <chromosome col>
#       tss: chrom_gene                           # (empty = the 'chrom' column;
#                                                 #  input names, as in the file)

LIFTOVER = config.get('liftover', {})
LIFTOVER_FROM = LIFTOVER.get('from_build', VARIANTS_BUILD)
LIFTOVER_TO = LIFTOVER.get('to_build', 'hg38')
LIFTOVER_CHAIN = LIFTOVER.get('chain')

LIFTOVER_COLUMN_ARGS = column_args(LIFTOVER.get('columns'), 'liftover.columns')
LIFTOVER_EXTRA_ARGS = " ".join(
    "--extra-pos-col " + shlex.quote(f"{pos}:{chrom}" if chrom else str(pos))
    for pos, chrom in (LIFTOVER.get('extra_pos_columns') or {}).items())

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
            chain=f"--chain {LIFTOVER_CHAIN}" if LIFTOVER_CHAIN else "",
            columns=LIFTOVER_COLUMN_ARGS,
            extra=LIFTOVER_EXTRA_ARGS,
        shell:
            """
            python {params.script} {input.variants} {output.lifted} \
                --from-build {params.from_build} --to-build {params.to_build} \
                 {params.chain} {params.columns} {params.extra}
            """

else:
    VARIANTS_HG38 = VARIANTS   # already in the target build, nothing to lift