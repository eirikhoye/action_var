# Download ensembl-vep:
cd resources/
singularity pull ensembl-vep.sif docker://ensemblorg/ensembl-vep:release_111.0

# Test if it works:
singularity exec ensembl-vep.sif vep --help

VEP description:
Ensembl VEP (Variant Effect Predictor) is a tool that takes genomic variants — SNVs, indels, structural variants — and annotates them with their predicted functional consequences.

Core function

Given a variant (e.g. from a VCF file), VEP determines which genes and transcripts it overlaps and predicts what effect it has, such as:

Missense, nonsense, synonymous, frameshift
Splice site disruption (donor/acceptor, splice region)
5'/3' UTR variants
Intronic, intergenic
Regulatory region and non-coding RNA effects
Upstream/downstream gene variants

For coding variants, it gives HGVS notation (c., p.), amino acid changes, and codon changes.

Consequence ranking

VEP assigns a "most severe consequence" per variant using Sequence Ontology terms, and can report all overlapping transcript consequences, not just the worst one — relevant when a variant hits multiple transcript isoforms differently.

Annotation sources it can pull in






# Modules

## Variant liftover and VEP annotation

Two scripts in `workflow/scripts/` take a variant table in hg19 coordinates, lift it to hg38, and annotate it with Ensembl VEP.

```
variants (hg19) ──► genome_liftover.py ──► variants (hg38) ──► annotate_vep.py ──► VEP-annotated VCF
                           │                                          │
                           └──► <output>.unmapped                     └──► <out>.input.vcf, <out>_summary.html
```

### Input format

Both scripts read a delimited text file (plain or `.gz`) with at least these columns:

```
VariantID;Chromosome;Position;Ref;Alt
10_89624151_C_T;10;89624151;C;T
```

- `Position` is 1-based (VCF convention).
- `Chromosome` may be written with or without a `chr` prefix.
- The separator (`;`, tab, `,` or `|`) is detected from the header line. Use `--sep` to set it explicitly, and quote it in the shell: `--sep ';'`, `--sep '\t'`.
- Extra columns are allowed and are kept by `genome_liftover.py`.

### `genome_liftover.py`: lift coordinates between builds

Converts `Position` from one genome build to another with UCSC chain files via [pyliftover](https://github.com/konstantint/pyliftover).

```bash
python workflow/scripts/genome_liftover.py \
    data/new_bridges_data/bridges_34genes.csv \
    data/new_bridges_data/bridges_34genes_hg38.csv
```

| Option | Default | Description |
|---|---|---|
| `--from-build` | `hg19` | Source build |
| `--to-build` | `hg38` | Target build |
| `--chain` | none | Local chain file, e.g. `hg19ToHg38.over.chain.gz`. Overrides the build options and avoids downloading from UCSC. |
| `--sep` | detected | Column separator |

**Output**

- `<output>`: same columns, order and separator as the input. `Position` holds the new coordinate, and `VariantID` is rebuilt as `<Chromosome>_<Position>_<Ref>_<Alt>`.
- `<output>.unmapped`: variants that could not be lifted, in their original coordinates. This file is only written if at least one variant fails.

A variant is treated as unmapped if:
1. it lies in a region with no chain alignment,
2. it lifts to a different chromosome,
3. it lifts to the minus strand (Ref/Alt are not reverse-complemented), or
4. its `Position` is missing.

**Caveats**

- Only the start position is lifted. This is exact for SNVs, but indels spanning a chain gap are not detected.
- Ref/Alt are not checked against the hg38 reference.
- Without `--chain`, the chain file is downloaded and cached in `~/.pyliftover`. Pass `--chain` to pin the exact file used.

### `annotate_vep.py`: annotate with VEP

Writes the variants to a minimal sites-only VCF (VCFv4.2) and runs Ensembl VEP 111 offline in a Singularity container against the local RefSeq cache. **The input must be in GRCh38 coordinates**, which is the output of `genome_liftover.py`.

```bash
python workflow/scripts/annotate_vep.py \
    data/new_bridges_data/bridges_34genes_hg38.csv \
    results/bridges_34genes_vep.vcf
```

| Option | Default | Description |
|---|---|---|
| `--root` | `/storage/mathelierarea/processed/eirikhoy/action_var` | Project root |
| `--vep-img` | `<root>/containers/ensembl-vep.sif` | VEP Singularity image |
| `--vep-cache` | `<root>/vep_cache` | Folder containing `homo_sapiens_refseq/` |
| `--sep` | detected | Column separator |

**Output**

- `<vep_out>`: VEP output in VCF format, with annotations in the `CSQ` INFO field.
- `<vep_out stem>.input.vcf`: the intermediate VCF passed to VEP, kept for provenance.
- `<vep_out>_summary.html`: VEP's run summary.

Rows with a missing `Position`, `Ref` or `Alt` are dropped, and the count is reported on stderr. The output folder is created if needed.

**VEP settings:** `--cache --offline --refseq --use_given_ref --assembly GRCh38 --vcf --force_overwrite`

- VEP looks for the cache at `<vep_cache>/homo_sapiens_refseq/111_GRCh38/`. The folder names must not be changed, and the cache version must match the VEP version in the container.
- The RefSeq cache is BAM-edited, and no FASTA is given, so `--use_given_ref` is required. The REF alleles are therefore **not** validated against the genome, and HGVS notation is unavailable. Add `--fasta` with a GRCh38 FASTA if either is needed.
- The container is started with `-B /storage`, so all input, output, image and cache paths must lie under `/storage`.

### Required resources

```
<root>/
├── containers/ensembl-vep.sif        # Ensembl VEP 111
└── vep_cache/
    └── homo_sapiens_refseq/
        ├── 111_GRCh37/
        └── 111_GRCh38/
```

Python dependencies: `pandas`, `numpy`, `pyliftover`, and `pytest` for the tests.

### Tests

```bash
pytest tests/ -v
```

The tests in `tests/test_scripts.py` need no network, Singularity or VEP cache. Liftover is tested against a small synthetic chain file with known outcomes: an offset, a chain gap, a chromosome change and a minus-strand hit. VEP is replaced by a fake `singularity` executable that records the command line, so the exact VEP invocation and the generated VCF can be checked.