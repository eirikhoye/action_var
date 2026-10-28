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

## Input variants

The `variants` entry in the config file points to a delimited text file (plain or `.gz`) with a header line and one variant per row. These columns are required:

| Column | Description |
|---|---|
| `VariantID` | Variant identifier, conventionally `<Chromosome>_<Position>_<Ref>_<Alt>` |
| `Chromosome` | Chromosome, with or without a `chr` prefix (`10` or `chr10`) |
| `Position` | 1-based position (VCF convention) |
| `Ref` | Reference allele (A/C/G/T) |
| `Alt` | Alternate allele (A/C/G/T) |

Example:

```
VariantID;Chromosome;Position;Ref;Alt
10_89624151_C_T;10;89624151;C;T
10_89624278_GA_G;10;89624278;GA;G
```

- **Separator:** `;`, tab, `,` or `|`, detected from the header line.
- **Genome build:** set `variants_build` in the config (e.g. `hg19`). If it is not `hg38`, the variants are lifted to hg38 before annotation. Variants that cannot be lifted are written to `<output>.unmapped`.
- **Indels:** use VCF-style representation, where Ref and Alt share the first base (e.g. `GA>G`, `A>AC`). Only the start position is lifted, so long indels are less reliable than SNVs.
- **Extra columns** are allowed and are passed through unchanged.
- **Recommended:** remove duplicate rows and left-normalise the variants (e.g. with `bcftools norm`) beforehand. Ref alleles are not checked against the reference genome.
