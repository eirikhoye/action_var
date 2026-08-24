from io import StringIO

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pandas as pd
from tqdm import tqdm
import argparse


API_KEY = "***REMOVED***"
dna_model = dna_client.create(API_KEY)


def run_alphagenome(args):
    vcf = pd.read_csv(args.input, sep='\t')

    required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']

    for column in required_columns:
        if column not in vcf.columns:
            raise ValueError(f'VCF file is missing required column: {column}.')

    organism = args.organism
    sequence_length = args.sequence_length
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]

    # Scorers
    score_rna_seq = args.score_rna_seq
    score_cage = args.score_cage
    score_procap = args.score_procap
    score_atac = args.score_atac
    score_dnase = args.score_dnase
    score_chip_histone = args.score_chip_histone
    score_chip_tf = args.score_chip_tf
    score_polyadenylation = args.score_polyadenylation
    score_splice_sites = args.score_splice_sites
    score_splice_site_usage = args.score_splice_site_usage
    score_splice_junctions = args.score_splice_junctions

    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS
    }
    organism = organism_map[organism]

    # Parse scorer specification.
    scorer_selections = {
        'rna_seq': score_rna_seq,
        'cage': score_cage,
        'procap': score_procap,
        'atac': score_atac,
        'dnase': score_dnase,
        'chip_histone': score_chip_histone,
        'chip_tf': score_chip_tf,
        'polyadenylation': score_polyadenylation,
        'splice_sites': score_splice_sites,
        'splice_site_usage': score_splice_site_usage,
        'splice_junctions': score_splice_junctions,
    }

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]

    # Drop scorers or outputs not supported for selected organism
    unsupported_scorers = [
        scorer
        for scorer in selected_scorers
            if (
                organism.value
                not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
            )
            | (
                (scorer.requested_output == dna_client.OutputType.PROCAP)
                & (organism == dna_client.Organism.MUS_MUSCULUS)
            )
    ]
    if len(unsupported_scorers) > 0:
        print(
            f'Excluding {unsupported_scorers} scorers as they are not supported for'
            f' {organism}.'
        )
        for unsupported_scorer in unsupported_scorers:
            selected_scorers.remove(unsupported_scorers)

    # Score variants in the VCF file
    results = []

    for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):
        variant = genome.Variant(
            chromosome=str(vcf_row.CHROM),
            position=int(vcf_row.POS),
            reference_bases=vcf_row.REF,
            alternate_bases=vcf_row.ALT,
            name=vcf_row.variant_id,
        )
        interval = variant.reference_interval.resize(sequence_length)

        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism,
        )
        results.append(variant_scores)

    df_scores = variant_scorers.tidy_scores(results)

    df_scores.to_csv(args.output, index=False)


def main():
    parser = argparse.ArgumentParser(
        description='predict functional score on variants with alphagenome'
    )

    # Required posiitonal arguments
    parser.add_argument(
        "input", 
        help='input tsv in VCF format for AlphaGenome'
    )
    parser.add_argument(
        'output', 
        help='output matrix'
    )

    # General settings
    parser.add_argument(
        "--organism", 
        default='human', 
        choices=['human', 'mouse'],
        help='human or mouse genome, (default: %(default)s)'
    )

    parser.add_argument(
        "--sequence_length", 
        default="1MB", 
        choices=['16KB', '100KB', '500KB', '1MB'],
        help="sequence around variants to predict, (default: %(default)s)")

    # ------------------------------------------------------------------
    # Scorers: ON by default, --no-X disables them
    # ------------------------------------------------------------------

    parser.add_argument(
        "--no-rna-seq",
        dest="score_rna_seq",
        action="store_false",
        help="Disable RNA-seq scorer",
    )

    parser.add_argument(
        "--no-cage",
        dest="score_cage",
        action="store_false",
        help="Disable CAGE scorer",
    )

    parser.add_argument(
        "--no-procap",
        dest="score_procap",
        action="store_false",
        help="Disable PRO-cap scorer",
    )

    parser.add_argument(
        "--no-atac",
        dest="score_atac",
        action="store_false",
        help="Disable ATAC scorer",
    )

    parser.add_argument(
        "--no-dnase",
        dest="score_dnase",
        action="store_false",
        help="Disable DNase scorer",
    )

    parser.add_argument(
        "--no-chip-histone",
        dest="score_chip_histone",
        action="store_false",
        help="Disable ChIP histone scorer",
    )

    parser.add_argument(
        "--no-chip-tf",
        dest="score_chip_tf",
        action="store_false",
        help="Disable ChIP TF scorer",
    )

    parser.add_argument(
        "--no-polyadenylation",
        dest="score_polyadenylation",
        action="store_false",
        help="Disable polyadenylation scorer",
    )

    parser.add_argument(
        "--no-splice-sites",
        dest="score_splice_sites",
        action="store_false",
        help="Disable splice sites scorer",
    )

    parser.add_argument(
        "--no-splice-site-usage",
        dest="score_splice_site_usage",
        action="store_false",
        help="Disable splice site usage scorer",
    )

    parser.add_argument(
        "--no-splice-junctions",
        dest="score_splice_junctions",
        action="store_false",
        help="Disable splice junctions scorer",
    )

    # All scorers default to True
    parser.set_defaults(
        score_rna_seq=True,
        score_cage=True,
        score_procap=True,
        score_atac=True,
        score_dnase=True,
        score_chip_histone=True,
        score_chip_tf=True,
        score_polyadenylation=True,
        score_splice_sites=True,
        score_splice_site_usage=True,
        score_splice_junctions=True,
    )


    args = parser.parse_args()

    run_alphagenome(args)


if __name__ == "__main__":
    main()