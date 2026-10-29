"""Shared helpers for the AlphaGenome scripts (API key, variant tables, retries).

The API key is never stored in the repository. It is read, in order, from
``--api-key-file`` (if given), the ``ALPHAGENOME_API_KEY`` environment variable,
or ``<project root>/resources/API_KEY`` (git-ignored, as in the notebooks).
"""

import os
import sys
import time
from pathlib import Path

import pandas as pd

DEFAULT_KEY_FILE = Path('resources/API_KEY')
VARIANT_COLUMNS = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']


def load_api_key(key_file=None):
    if key_file:
        return Path(key_file).read_text().strip()
    if os.environ.get('ALPHAGENOME_API_KEY'):
        return os.environ['ALPHAGENOME_API_KEY'].strip()
    if DEFAULT_KEY_FILE.exists():
        return DEFAULT_KEY_FILE.read_text().strip()
    sys.exit("No AlphaGenome API key: set ALPHAGENOME_API_KEY, pass --api-key-file, "
             f"or put it in {DEFAULT_KEY_FILE} (git-ignored)")


def read_variants(paths):
    """Concatenate variant TSVs (variant_id, CHROM, POS, REF, ALT [, gene]); skips empty files."""
    frames = []
    for path in paths:
        df = pd.read_csv(path, sep='\t', dtype={'CHROM': str})
        missing = [c for c in VARIANT_COLUMNS if c not in df]
        if missing:
            raise ValueError(f"{path} is missing column(s) {missing}")
        frames.append(df)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=VARIANT_COLUMNS)
    df['CHROM'] = 'chr' + df.CHROM.astype(str).str.removeprefix('chr')
    df['POS'] = df.POS.astype(int)
    return df.drop_duplicates('variant_id').reset_index(drop=True)


def with_retries(fn, what, attempts=5, wait=10):
    """Call ``fn()``; retry transient (gRPC/network) errors with linear back-off.

    ValueError (bad request / not found) and PermissionError (bad API key) are
    raised at once: retrying cannot fix them.
    """
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except (ValueError, PermissionError):
            raise
        except Exception as e:  # noqa: BLE001 - the client raises assorted gRPC errors
            if attempt == attempts:
                raise RuntimeError(f"{what} failed after {attempts} attempts") from e
            print(f"retry {attempt}/{attempts} {what}: {type(e).__name__}: {str(e)[:200]}",
                  file=sys.stderr)
            time.sleep(wait * attempt)


def chunks(df, size):
    for start in range(0, len(df), size):
        yield start // size, df.iloc[start:start + size]
