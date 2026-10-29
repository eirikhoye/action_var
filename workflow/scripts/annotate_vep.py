#!/usr/bin/env python3
"""
annotate_vep.py
===============
 
Annotate a delimited variant table with Ensembl VEP (run inside a
Singularity container, offline, against a local RefSeq cache).
 
Purpose
-------
Converts a variant table into a minimal sites-only VCF and runs VEP on it,
producing a VEP-annotated VCF. Intended to run *after* genome_liftover.py,
since VEP is called with ``--assembly GRCh38``.
 
Input
-----
A delimited text file (optionally gzipped) with at least these columns::
 
    VariantID;Chromosome;Position;Ref;Alt

* Other column names can be given with ``--id-col``, ``--chrom-col``,
  ``--pos-col``, ``--ref-col`` and ``--alt-col``.
* Coordinates must be **GRCh38**, 1-based.
* ``Chromosome`` may be given with or without a ``chr`` prefix.
* Rows with a missing ``Position``, ``Ref`` or ``Alt`` are dropped, and the
  number dropped is reported on stderr.
* The separator is detected from the header line unless ``--sep`` is given.
 
Output
------
``<vep_out>``
    VEP output in VCF format; annotations are in the ``CSQ`` INFO field.
``<vep_out stem>.input.vcf``
    The intermediate VCFv4.2 file passed to VEP (kept for provenance).
``<vep_out>_summary.html``
    VEP's run summary, written by VEP itself.
 
The output directory is created if it does not exist.
 
VEP configuration
-----------------
VEP is run with::
 
    --cache --offline --refseq --use_given_ref --assembly GRCh38
    --mane --format vcf --vcf --force_overwrite
 
* ``--refseq`` selects the cache in ``<vep_cache>/homo_sapiens_refseq/``.
  VEP resolves the full path as
  ``<dir_cache>/homo_sapiens_refseq/<VEP version>_GRCh38/``, so the cache
  version must match the VEP version in the container (here 111).
* ``--use_given_ref`` is required because the RefSeq cache is BAM-edited and
  no FASTA is supplied. VEP trusts the REF alleles in the input; they are
  **not** validated against the genome, and HGVS notation is unavailable.
  Supply a GRCh38 FASTA (``--fasta``) instead if either matters.
* ``--mane`` adds the ``MANE_SELECT`` and ``MANE_PLUS_CLINICAL`` CSQ
  fields, which vep_report.py uses to pick one transcript per gene.
* VCF records are written sorted by chromosome name (lexicographic) and
  position. VEP does not require sorted input; the order only aids reading.
* ``ID`` is set to ``VariantID`` (``.`` if missing); QUAL, FILTER and INFO
  are ``.``.
 
Plugins (optional)
------------------
``--sift``/``--polyphen`` read predictions from the cache. ``--plugin SPEC``
(repeatable) is passed to VEP unchanged; the plugin *code* ships with the
ensembl-vep image, but plugin *data* does not and must be on a bound path.
The ``annotate_vep`` rule builds these from ``vep.plugins`` in the config;
the CSQ fields they add are read by classify_variant_effects.py:

=============  ===========================================  =======================
Plugin         Spec                                         CSQ fields
=============  ===========================================  =======================
dbNSFP         ``dbNSFP,<file>,REVEL_score,MetaRNN_score,``  same as the field list
               ``AlphaMissense_score,AlphaMissense_pred``
AlphaMissense  ``AlphaMissense,file=<tsv.gz>``              am_class,
                                                            am_pathogenicity
REVEL          ``REVEL,file=<tsv.gz>``                      REVEL
dbscSNV        ``dbscSNV,<file>``                           ada_score, rf_score
MaxEntScan     ``MaxEntScan,<fordownload dir>``             MaxEntScan_ref/alt/diff
SpliceAI       ``SpliceAI,snv=<vcf.gz>,indel=<vcf.gz>``     SpliceAI_pred_DS_AG/AL/
                                                            DG/DL
=============  ===========================================  =======================

Containers and paths
--------------------
The container is run with ``-B /storage``. Input, output, image and cache
paths must therefore resolve under ``/storage`` (or another path already
visible in the container, e.g. the current working directory).
 
Usage
-----
::
 
    python annotate_vep.py INPUT VEP_OUT [--sep SEP] [--root DIR]
                           [--vep-img SIF] [--vep-cache DIR]
                           [--id-col COL] [--chrom-col COL] [--pos-col COL]
                           [--ref-col COL] [--alt-col COL]
 
Defaults: ``--vep-img <root>/containers/ensembl-vep.sif``,
``--vep-cache <root>/vep_cache``.
 
Example
-------
::
 
    python workflow/scripts/annotate_vep.py \\
        test/bridges_34genes_hg38.csv test/bridges_34genes_vep.vcf
 
Dependencies
------------
Python >= 3.9, pandas, Singularity/Apptainer, an Ensembl VEP image and a
matching offline RefSeq cache.
 
Author: Eirik Høye
"""

import argparse
import codecs
import gzip
import subprocess
import sys
from pathlib import Path
 
import pandas as pd

REQUIRED_COLUMNS = ["VariantID", "Chromosome", "Position", "Ref", "Alt"]


def parse_args(argv=None):
    """Parse command-line arguments and fill path defaults from ``--root``."""
    parser = argparse.ArgumentParser(
        description="Annotate variants with VEP",
        epilog="See the module docstring for VEP settings and assumptions.",
    )
    parser.add_argument("variants_csv",
                        help="Input file with VariantID;Chromosome;Position;Ref;Alt (GRCh38)")
    parser.add_argument("vep_out", type=Path, help="VEP output file (VCF)")
    parser.add_argument("--sep", default=None,
                        help=r"Column separator, e.g. ';' or '\t' (default: detected from header)")
    parser.add_argument("--root", type=Path,
                        default=Path("/storage/mathelierarea/processed/eirikhoy/action_var"),
                        help="Project root dir")
    parser.add_argument("--vep-img", type=Path, default=None,
                        help="VEP singularity image (default: <root>/containers/ensembl-vep.sif)")
    parser.add_argument("--vep-cache", type=Path, default=None,
                        help="VEP cache dir containing homo_sapiens_refseq/ "
                             "(default: <root>/vep_cache)")
    parser.add_argument("--id-col",    default="VariantID",  help="Variant ID column (default: %(default)s)")
    parser.add_argument("--chrom-col", default="Chromosome", help="Chromosome column (default: %(default)s)")
    parser.add_argument("--pos-col",   default="Position",   help="Position column (default: %(default)s)")
    parser.add_argument("--ref-col",   default="Ref",        help="Ref allele column (default: %(default)s)")
    parser.add_argument("--alt-col",   default="Alt",        help="Alt allele column (default: %(default)s)")

    # Optional effect-prediction annotations (see "Plugins" in the module docstring)
    parser.add_argument("--sift", choices=["p", "s", "b"], default=None,
                        help="Add SIFT from the cache: p=prediction, s=score, b=both")
    parser.add_argument("--polyphen", choices=["p", "s", "b"], default=None,
                        help="Add PolyPhen from the cache: p=prediction, s=score, b=both")
    parser.add_argument("--plugin", action="append", default=[], metavar="SPEC",
                        help="VEP plugin spec, passed as --plugin SPEC (repeatable), "
                             "e.g. 'dbscSNV,/storage/.../dbscSNV1.1_GRCh38.txt.gz'")
    parser.add_argument("--dir-plugins", type=Path, default=None,
                        help="Plugin code directory inside the container (default: VEP's own)")
    parser.add_argument("--bind", action="append", default=[], metavar="PATH",
                        help="Extra Singularity bind path, e.g. for plugin data outside /storage "
                             "(repeatable)")

    args = parser.parse_args(argv)

    # Fill path defaults from root
    args.vep_img = args.vep_img or args.root / "containers/ensembl-vep.sif"
    args.vep_cache = args.vep_cache or args.root / "vep_cache"
    return args


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


def prepare_variants(variants):
    """Clean a variant table for VCF output"""
    variants = variants.copy()
    variants["chrom"] = "chr" + variants["Chromosome"].astype(str).str.removeprefix("chr")
    n_before = len(variants)
    variants = variants.dropna(subset=["Position", "Ref", "Alt"])
    n_dropped = n_before - len(variants)
    if n_dropped:
        print(f"Dropped {n_dropped}/{n_before} variants with missing "
              "Position, Ref or Alt", file=sys.stderr)
    return (variants.astype({"Position": int})
                    .sort_values(["chrom", "Position"]))


def write_vcf(variants, vcf_path):
    """Write a minimal sites-only VCFv4.2 file."""    
    with open(vcf_path, 'w') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for v in variants.itertuples():
            vid = "." if pd.isna(v.VariantID) else v.VariantID
            f.write(f"{v.chrom}\t{v.Position}\t{vid}\t{v.Ref}\t{v.Alt}\t.\t.\t.\n")

            
def vep_extra_args(sift=None, polyphen=None, plugins=(), dir_plugins=None):
    """VEP arguments for the optional SIFT/PolyPhen and plugin annotations."""
    extra = []
    if sift:
        extra += ['--sift', sift]
    if polyphen:
        extra += ['--polyphen', polyphen]
    if dir_plugins:
        extra += ['--dir_plugins', str(dir_plugins)]
    for spec in plugins:
        extra += ['--plugin', spec]
    return extra


def annotate_vep(vcf_in, vep_out, vep_img, vep_cache, extra_args=(), binds=()):
    """Run VEP in Singularity on ``vcf_in`` and write ``vep_out``.

    Raises
    ------
    subprocess.CalledProcessError
        If VEP exits with a non-zero status.
    """
    bind_args = []
    for path in ['/storage', *binds]:
        bind_args += ['-B', str(path)]
    cmd = [
        'singularity', 'exec', *bind_args, str(vep_img), 'vep',
        '--input_file', str(vcf_in),
        '--format', 'vcf',
        '--output_file', str(vep_out),
        '--vcf',
        '--cache', '--offline', '--refseq', '--use_given_ref',
        '--dir_cache', str(vep_cache),
        '--assembly', 'GRCh38',
        '--mane',
        '--force_overwrite',
        *extra_args,
    ]
    print('Running:', ' '.join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


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

    # Only the canonical columns are needed for the VCF; this also avoids
    # clashes with unmapped input columns (e.g. one already named "chrom")
    variants = pd.read_csv(args.variants_csv, sep=sep).rename(columns=column_map)
    check_columns(variants, args.variants_csv, column_map)
    variants = prepare_variants(variants[REQUIRED_COLUMNS])

    args.vep_out.parent.mkdir(parents=True, exist_ok=True)
    vcf_in = args.vep_out.with_suffix(".input.vcf")
    write_vcf(variants, vcf_in)
    extra = vep_extra_args(args.sift, args.polyphen, args.plugin, args.dir_plugins)
    annotate_vep(vcf_in, args.vep_out, args.vep_img, args.vep_cache, extra, args.bind)


if __name__ == "__main__":
    main()