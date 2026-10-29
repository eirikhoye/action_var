# rules/alphagenome_ism.smk
# DESIGN SKELETON: in silico mutagenesis (ISM) of the non-coding genome around the
# 34 BRIDGES genes with the AlphaGenome Atlas, then Sei on the distal elements that
# AlphaGenome predicts to affect the gene. See docs/ism_sei_design.md.
#
# Status:
#   ism_regions, ism_atlas          runnable (Atlas interval queries, cached, per gene)
#   ism_elements, ism_sei_*, report stubs: raise NotImplementedError when run
#
# Requires annotate_vep (for vep.report.genes / gtf). All keys are optional except where noted:
#
#   alphagenome_ism:
#     outdir: results/ism
#     genes: [BRCA1, BRCA2]         # default: vep.report.genes
#     flank: 100000                 # bp each side of the MANE transcript
#     around_tss: false             # true: window = TSS +/- flank instead
#     chunk: 4000                   # bp per Atlas query (one cached parquet each)
#     scorers: [AVI_SCORE, RNA_SEQ, SPLICE_SITES, DNASE, CAGE, CHIP_HISTONE]
#     ontology_terms: [CL:0002327]
#     elements:                     # ism_elements (to implement)
#       score: "RNA_SEQ:CL:0002327 polyA plus RNA-seq"   # gene-level effect track
#       min_quantile: 0.99          # Atlas calibrated quantile of |effect|
#       exclude: [CDS]              # region classes to drop before calling elements
#       merge_distance: 50          # bp between hit positions merged into one element
#       min_positions: 3
#     sei:
#       top_snvs_per_element: 10    # strongest ISM SNVs per element sent to Sei VEP
#       seqclass_annotation: ...    # hg38 Sei sequence-class BED (Zenodo 7113988)

if 'VEP_REPORT' not in globals():
    raise WorkflowError("alphagenome_ism needs the annotate_vep module (vep.report.genes/gtf)")

ISM = config.get('alphagenome_ism') or {}
ISM_DIR = ISM.get('outdir', 'results/ism')
ISM_GENES = ISM.get('genes') or VEP_REPORT['genes']
ISM_REGIONS = f"{ISM_DIR}/regions.tsv"
_ISM_ELEM = ISM.get('elements') or {}
_ISM_SEI = ISM.get('sei') or {}
_ISM_KEY = (config.get('alphagenome') or {}).get('api_key_file')

ALL_OUTPUTS.append(ISM_REGIONS)
ALL_OUTPUTS.extend(expand(f"{ISM_DIR}/atlas/{{gene}}.parquet", gene=ISM_GENES))
if ISM.get('run_design_steps'):   # stubs below; switch on once implemented
    ALL_OUTPUTS.append(f"{ISM_DIR}/ism_sei_report.html")


# 1. One window per gene around its MANE transcript.
rule ism_regions:
    output:
        regions=ISM_REGIONS
    params:
        script="workflow/scripts/make_ism_regions.py",
        gtf=VEP_REPORT['gtf'],
        genes=','.join(ISM_GENES),
        flank=ISM.get('flank', 100_000),
        tss="--around-tss" if ISM.get('around_tss') else "",
    shell:
        """
        python {params.script} {params.gtf} {output.regions} --genes {params.genes} \
            --flank {params.flank} {params.tss}
        """


# 2. Every SNV in the window from the Atlas, effects on that gene (wide table).
rule ism_atlas:
    input:
        regions=ISM_REGIONS
    output:
        scores=f"{ISM_DIR}/atlas/{{gene}}.parquet"
    params:
        script="workflow/scripts/query_alphagenome_atlas.py",
        scorers=','.join(ISM.get('scorers', ['AVI_SCORE', 'RNA_SEQ', 'SPLICE_SITES', 'DNASE',
                                            'CAGE', 'CHIP_HISTONE'])),
        ontology=shlex.quote(','.join(ISM.get('ontology_terms', ['CL:0002327']))),
        chunk=ISM.get('chunk', 4000),
        cache=f"{ISM_DIR}/atlas/.cache",
        key=f"--api-key-file {shlex.quote(_ISM_KEY)}" if _ISM_KEY else "",
    shell:
        """
        python {params.script} interval {input.regions} {output.scores} --gene {wildcards.gene} \
            --scorers {params.scorers} --ontology-terms {params.ontology} \
            --chunk {params.chunk} --cache {params.cache} {params.key}
        """


# 3. (to implement) Call distal elements: per position the max-|effect| ALT on the gene's
#    RNA_SEQ track; keep positions with quantile >= min_quantile outside excluded regions;
#    merge hits within merge_distance; keep elements with >= min_positions.
#    Output BED: chrom, start, end, gene, n_positions, max_effect, direction, tss_distance.
rule ism_elements:
    input:
        scores=f"{ISM_DIR}/atlas/{{gene}}.parquet",
        regions=ISM_REGIONS,
    output:
        bed=f"{ISM_DIR}/elements/{{gene}}.bed",
        top_snvs=f"{ISM_DIR}/elements/{{gene}}_top_snvs.tsv",
    run:
        raise NotImplementedError("ism_elements: design only, see docs/ism_sei_design.md")


# 4a. (to implement) What each element *is*: overlap with Sei's hg38 sequence-class
#     annotation (Zenodo 7113988). No model run needed.
rule ism_sei_annotate_elements:
    input:
        bed=f"{ISM_DIR}/elements/{{gene}}.bed"
    output:
        tsv=f"{ISM_DIR}/sei/{{gene}}_element_classes.tsv"
    run:
        raise NotImplementedError("ism_sei_annotate_elements: design only")


# 4b. (to implement) What the strongest ISM mutations *do* to chromatin: the top SNVs per
#     element through the same Sei steps as rules/sei.smk (prepare VCF -> sei_chromatin ->
#     sei_seqclass -> merge), so AlphaGenome RNA effect and Sei class change can be compared.
rule ism_sei_variants:
    input:
        snvs=f"{ISM_DIR}/elements/{{gene}}_top_snvs.tsv"
    output:
        tsv=f"{ISM_DIR}/sei/{{gene}}_top_snvs_sei.tsv"
    run:
        raise NotImplementedError("ism_sei_variants: design only")


# 5. (to implement) Per gene: element map, AlphaGenome vs Sei concordance, and overlap
#    with the observed BRIDGES variants from variant_effects.
rule ism_sei_report:
    input:
        classes=expand(f"{ISM_DIR}/sei/{{gene}}_element_classes.tsv", gene=ISM_GENES),
        snvs=expand(f"{ISM_DIR}/sei/{{gene}}_top_snvs_sei.tsv", gene=ISM_GENES),
    output:
        html=f"{ISM_DIR}/ism_sei_report.html"
    run:
        raise NotImplementedError("ism_sei_report: design only")
