#!/usr/bin/env python3
"""
genome_liftover.py
==================
 
Lift single-nucleotide variant coordinates between genome builds
(default: hg19/GRCh37 -> hg38/GRCh38) using UCSC chain files via pyliftover.
 
Purpose
-------
Converts a delimited variant table so that downstream tools (e.g. VEP run
with ``--assembly GRCh38``) receive coordinates in the correct build. The
output has the same columns, column order and separator as the input; only
``Position`` and ``VariantID`` are rewritten.
 
Input
-----
A delimited text file (optionally gzipped) with at least these columns::
 
    VariantID;Chromosome;Position;Ref;Alt
 
* ``Chromosome`` may be given with or without a ``chr`` prefix (``10`` or
  ``chr10``). The prefix style of the input is preserved in the output.
* ``Position`` is 1-based (VCF convention).
* The separator is detected from the header line unless ``--sep`` is given.
* Other column names can be given with ``--id-col``, ``--chrom-col``,
  ``--pos-col``, ``--ref-col`` and ``--alt-col``. The output keeps the input's
  column names.
* Further coordinate columns (e.g. an end position or a gene TSS) can be
  lifted with ``--extra-pos-col COL[:CHROM_COL]`` (repeatable). ``CHROM_COL``
  defaults to the chromosome column.

Output
------
``<output>``
    Successfully lifted variants. ``Position`` holds the target-build
    coordinate and ``VariantID`` is rebuilt as ``<Chromosome>_<Position>_<Ref>_<Alt>``.
``<output>.unmapped``
    Variants that could not be lifted, in their *original* coordinates.
    Only written if at least one variant fails.
 
A one-line summary (``N/M variants failed to lift``) is printed to stderr.
 
Method and coordinate conventions
---------------------------------
* pyliftover works in 0-based coordinates; positions are converted
  1-based -> 0-based before lifting and back to 1-based afterwards.
* Only the first hit returned by pyliftover is used.
* A variant is treated as *unmapped* (and written to ``.unmapped``) if:
 
  1. it falls in a region with no chain alignment (deleted/rearranged
     between builds),
  2. it lifts to a different chromosome than the source, or
  3. it lifts to the minus strand. Ref/Alt would then have to be
     reverse-complemented, which this script deliberately does not do.
  4. its ``Position`` is missing.

* Extra position columns are lifted by the same rules, but a failure there
  does not drop the variant: the value is left empty, and the number of
  such values is reported on stderr.

* Only the start coordinate is lifted. This is exact for SNVs. For indels
  or MNVs spanning a chain gap the end coordinate may not be contiguous in
  the target build; such cases are not detected.
* Ref/Alt alleles are not checked against the target reference. Run VEP
  (or ``bcftools norm --check-ref``) with a FASTA to verify them.
 
Reproducibility
---------------
By default pyliftover downloads the chain file from UCSC and caches it in
``~/.pyliftover``. For offline or HPC use, and to pin the exact chain file
used, pass it explicitly with ``--chain`` (e.g. ``hg19ToHg38.over.chain.gz``).
 
Usage
-----
::
 
    python genome_liftover.py INPUT OUTPUT [--from-build hg19] [--to-build hg38]
                              [--chain FILE] [--sep SEP]
                              [--id-col COL] [--chrom-col COL] [--pos-col COL]
                              [--ref-col COL] [--alt-col COL]
                              [--extra-pos-col COL[:CHROM_COL] ...]

Examples
--------
::
 
    # hg19 -> hg38, separator detected from header
    python workflow/scripts/genome_liftover.py \\
        data/new_bridges_data/bridges_34genes.csv \\
        data/new_bridges_data/bridges_34genes_hg38.csv
 
    # Explicit separator (quote it, ';' ends a command in bash)
    python workflow/scripts/genome_liftover.py in.tsv out.tsv --sep '\\t'
 
    # Offline, with a local chain file
    python workflow/scripts/genome_liftover.py in.csv out.csv \\
        --chain resources/hg19ToHg38.over.chain.gz
 
Dependencies
------------
Python >= 3.9, pandas, numpy, pyliftover.
 
Author: Eirik Høye
"""

import argparse
import codecs
import gzip
import sys

import numpy as np
import pandas as pd
from pyliftover import LiftOver

REQUIRED_COLUMNS = ["VariantID", "Chromosome", "Position", "Ref", "Alt"]


def parse_args(argv=None):
    """Parse command-line arguments (``argv=None`` reads ``sys.argv``)."""
    parser = argparse.ArgumentParser(
        description="Lift variant coordinates between genome builds",
        epilog="See the module docstring for conventions and filtering rules.",
    )
    parser.add_argument("variants_csv",
                        help="Input file with VariantID;Chromosome;Position;Ref;Alt")
    parser.add_argument("output", help="Output file (same format as input)")
    parser.add_argument("--from-build", default="hg19",
                        help="Source genome build (default: hg19)")
    parser.add_argument("--to-build", default="hg38",
                        help="Target genome build (default: hg38)")
    parser.add_argument("--chain", default=None,
                        help="Local UCSC chain file. Overrides --from-build/--to-build "
                             "and avoids downloading from UCSC")
    parser.add_argument("--sep", default=None,
                        help=r"Column separator, e.g. ';' or '\t' "
                             "(default: detected from header line)")
    parser.add_argument("--id-col",    default="VariantID",  help="Variant ID column (default: %(default)s)")
    parser.add_argument("--chrom-col", default="Chromosome", help="Chromosome column (default: %(default)s)")
    parser.add_argument("--pos-col",   default="Position",   help="Position column (default: %(default)s)")
    parser.add_argument("--ref-col",   default="Ref",        help="Ref allele column (default: %(default)s)")
    parser.add_argument("--alt-col",   default="Alt",        help="Alt allele column (default: %(default)s)")
    parser.add_argument("--extra-pos-col", action="append", default=[],
                        metavar="COL[:CHROM_COL]",
                        help="Additional 1-based position column to lift, optionally with "
                             "its own chromosome column (default: --chrom-col). Repeatable")
    return parser.parse_args(argv)
    

def guess_sep(path, candidates=(";", "\t", ",", "|")):
    """Guess the column separator from the header line."""
    path = str(path)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.strip() and not line.startswith("##"):
                header = line
                break
        else:
            raise ValueError(f"{path} appears to be empty")

    counts = {sep: header.count(sep) for sep in candidates}
    best = max(counts, key=counts.get)
    if counts[best] == 0:
        raise ValueError(f"Could not detect separator in header of {path}: {header!r}")
    return best


def check_columns(df, path, column_map):
    """Raise ``ValueError`` if a required column is missing after renaming."""
    inverse = {v: k for k, v in column_map.items()}
    missing = [inverse[c] for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path} is missing column(s) {missing}; found {list(df.columns)}. "
            "Check the --*-col options and the separator.")


def lift_position(lift, chrom, pos):
    """Lift one 1-based position to the target build.
 
    Parameters
    ----------
    lift : pyliftover.LiftOver
        Initialised LiftOver object.
    chrom : str or int
        Chromosome, with or without ``chr`` prefix.
    pos : int or float
        1-based position in the source build. NaN is treated as unmapped.
 
    Returns
    -------
    int or float
        1-based position in the target build, or ``np.nan`` if the position
        is missing, unmappable, lifts to another chromosome, or lifts to the
        minus strand.
    """
    if pd.isna(pos):
        return np.nan
    chrom = "chr" + str(chrom).removeprefix("chr")
    hit = lift.convert_coordinate(chrom, int(pos) - 1)
    if not hit or hit[0][0] != chrom or hit[0][2] != "+":
        return np.nan
    return hit[0][1] + 1


def lift_column(lift, chroms, positions):
    """Lift a column of positions; returns a float Series with NaN for failures."""
    return pd.Series([lift_position(lift, c, p) for c, p in zip(chroms, positions)],
                     index=positions.index, dtype="float64")


def parse_extra_pos_cols(specs, column_map):
    """Parse ``COL[:CHROM_COL]`` specs into ``(pos_col, chrom_col)`` pairs.

    Names are given as in the input file and translated through
    ``column_map``, so the chromosome column may be referred to by its input name.
    """
    pairs = []
    for spec in specs:
        pos_col, _, chrom_col = spec.partition(":")
        if pos_col in column_map:
            raise ValueError(f"--extra-pos-col {pos_col!r} is already one of the --*-col columns")
        chrom_col = column_map.get(chrom_col, chrom_col) if chrom_col else "Chromosome"
        pairs.append((pos_col, chrom_col))
    return pairs


def main(argv=None):
    args = parse_args(argv)
    sep = codecs.decode(args.sep, "unicode_escape") if args.sep \
        else guess_sep(args.variants_csv)

    column_map = {
        args.id_col:    "VariantID",
        args.chrom_col: "Chromosome",
        args.pos_col:   "Position",
        args.ref_col:   "Ref",
        args.alt_col:   "Alt",
    }
    if len(column_map) != len(REQUIRED_COLUMNS):
        raise ValueError("The same input column was given for more than one --*-col option")
    inverse_map = {v: k for k, v in column_map.items()}

    variants = pd.read_csv(args.variants_csv, sep=sep).rename(columns=column_map)
    check_columns(variants, args.variants_csv, column_map)
    extra_pos_cols = parse_extra_pos_cols(args.extra_pos_col, column_map)
    missing = sorted({c for pair in extra_pos_cols for c in pair} - set(variants.columns))
    if missing:
        raise ValueError(f"--extra-pos-col: column(s) {missing} not found in {args.variants_csv}")

    lift = LiftOver(args.chain) if args.chain \
        else LiftOver(args.from_build, args.to_build)

    lifted = lift_column(lift, variants["Chromosome"], variants["Position"])
    failed = lifted.isna()

    print(f"{failed.sum()}/{len(variants)} variants failed to lift", file=sys.stderr)
    if failed.any():
        unmapped_path = args.output + ".unmapped"
        variants[failed].rename(columns=inverse_map).to_csv(unmapped_path, sep=sep, index=False)
        print(f"Unmapped variants written to {unmapped_path}", file=sys.stderr)

    variants = variants[~failed].copy()
    variants["Position"] = lifted[~failed].astype(int)
    variants['VariantID'] = (variants['Chromosome'].astype(str) + '_' 
                             + variants['Position'].astype(str) + '_' 
                             + variants['Ref'].astype(str) + '_' 
                             + variants['Alt'].astype(str))

    for pos_col, chrom_col in extra_pos_cols:
        lifted_extra = lift_column(lift, variants[chrom_col], variants[pos_col])
        n_failed = (lifted_extra.isna() & variants[pos_col].notna()).sum()
        print(f"{pos_col}: {n_failed}/{len(variants)} values failed to lift (left empty)",
              file=sys.stderr)
        variants[pos_col] = lifted_extra.astype("Int64")

    variants.rename(columns=inverse_map).to_csv(args.output, sep=sep, index=False)

if __name__ == "__main__":
    main()