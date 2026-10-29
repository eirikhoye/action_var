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


# Optional effect-prediction annotations (all keys optional; omit a plugin to skip it).
# Plugin data is not in the container: give paths under /storage, or add them to bind.
#
#   vep:
#     sift: b                      # p = prediction, s = score, b = both
#     polyphen: b
#     bind: [/cluster/projects/...] # extra Singularity bind paths
#     dir_plugins: /plugins        # only if VEP can't find the image's plugins
#     plugins:
#       dbNSFP:        {file: .../dbNSFP4.9a_grch38.gz,
#                       fields: [REVEL_score, MetaRNN_score, AlphaMissense_score, AlphaMissense_pred]}
#       AlphaMissense: {file: .../AlphaMissense_hg38.tsv.gz}
#       REVEL:         {file: .../new_tabbed_revel_grch38.tsv.gz}
#       dbscSNV:       {file: .../dbscSNV1.1_GRCh38.txt.gz}
#       MaxEntScan:    {dir: .../fordownload}
#       SpliceAI:      {snv: .../spliceai_scores.raw.snv.hg38.vcf.gz,
#                       indel: .../spliceai_scores.raw.indel.hg38.vcf.gz}
#     extra_plugins: []            # raw '--plugin' specs, passed through unchanged
_DBNSFP_DEFAULT_FIELDS = ['REVEL_score', 'MetaRNN_score', 'AlphaMissense_score', 'AlphaMissense_pred']


def vep_plugin_specs(vep):
    """(specs, data_files): VEP --plugin specs and the data paths they need."""
    plugins = vep.get('plugins') or {}
    known = {'dbNSFP', 'AlphaMissense', 'REVEL', 'dbscSNV', 'MaxEntScan', 'SpliceAI'}
    unknown = set(plugins) - known
    if unknown:
        raise ValueError(f"vep.plugins: unknown plugin(s) {sorted(unknown)}; "
                         f"allowed: {sorted(known)} (use vep.extra_plugins for others)")
    specs, files = [], []
    for name, opts in plugins.items():
        opts = opts or {}
        if name == 'dbNSFP':
            fields = opts.get('fields', _DBNSFP_DEFAULT_FIELDS)
            specs.append(','.join(['dbNSFP', opts['file'], *fields]))
            files.append(opts['file'])
        elif name in ('AlphaMissense', 'REVEL'):
            specs.append(f"{name},file={opts['file']}")
            files.append(opts['file'])
        elif name == 'dbscSNV':
            specs.append(f"dbscSNV,{opts['file']}")
            files.append(opts['file'])
        elif name == 'MaxEntScan':
            specs.append(f"MaxEntScan,{opts['dir']}")
            files.append(opts['dir'])
        elif name == 'SpliceAI':
            parts = [f"{k}={opts[k]}" for k in ('snv', 'indel') if opts.get(k)]
            if not parts:
                raise ValueError("vep.plugins.SpliceAI needs at least one of snv/indel")
            specs.append(','.join(['SpliceAI', *parts]))
            files.extend(opts[k] for k in ('snv', 'indel') if opts.get(k))
    specs.extend(vep.get('extra_plugins') or [])
    return specs, files


VEP_PLUGIN_SPECS, VEP_PLUGIN_FILES = vep_plugin_specs(VEP)
_vep_extra = []
for _key in ('sift', 'polyphen'):
    if VEP.get(_key):
        _vep_extra.append(f"--{_key} {VEP[_key]}")
if VEP.get('dir_plugins'):
    _vep_extra.append(f"--dir-plugins {shlex.quote(str(VEP['dir_plugins']))}")
_vep_extra += [f"--bind {shlex.quote(str(b))}" for b in VEP.get('bind') or []]
_vep_extra += [f"--plugin {shlex.quote(s)}" for s in VEP_PLUGIN_SPECS]
VEP_EXTRA_ARGS = " ".join(_vep_extra)


rule annotate_vep:
    input:
        var=VARIANTS_HG38,
        plugin_data=VEP_PLUGIN_FILES,   # fail early if plugin data is missing
    output:
        vcf=VARIANTS_VEP
    params:
        script="workflow/scripts/annotate_vep.py",
        columns=VEP_COLUMN_ARGS,
        extra=VEP_EXTRA_ARGS,
    shell:
        """
        python {params.script} {input.var} {output.vcf} {params.columns} {params.extra}
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
