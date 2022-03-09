"""
Function for reading the table created by the 'parse' subcommand.
"""
import logging

import numpy as np
import pandas as pd

from .utils import nt_to_aa

logger = logging.getLogger(__name__)


# This lists only the AIRR integer columns for now. The other columns, which are currently
# converted according to _STRING_COLUMNS or _INTEGER_COLUMNS should also be
# added here.
COLUMN_TYPES = {
    "v_alignment_start": "Int64",
    "v_alignment_end": "Int64",
    "d_alignment_start": "Int64",
    "d_alignment_end": "Int64",
    "j_alignment_start": "Int64",
    "j_alignment_end": "Int64",
    "v_sequence_start": "Int64",
    "v_sequence_end": "Int64",
    "v_germline_start": "Int64",
    "v_germline_end": "Int64",
    "d_sequence_start": "Int64",
    "d_sequence_end": "Int64",
    "d_germline_start": "Int64",
    "d_germline_end": "Int64",
    "j_sequence_start": "Int64",
    "j_sequence_end": "Int64",
    "j_germline_start": "Int64",
    "j_germline_end": "Int64",
    "fwr1_start": "Int64",
    "fwr1_end": "Int64",
    "cdr1_start": "Int64",
    "cdr1_end": "Int64",
    "fwr2_start": "Int64",
    "fwr2_end": "Int64",
    "cdr2_start": "Int64",
    "cdr2_end": "Int64",
    "fwr3_start": "Int64",
    "fwr3_end": "Int64",
    "fwr4_start": "Int64",
    "fwr4_end": "Int64",
    "cdr3_start": "Int64",
    "cdr3_end": "Int64",
}

# These columns contain string data
# convert them to str to avoid a PerformanceWarning
# TODO some of these are actually categorical or bool
_STRING_COLUMNS = [
    'v_call',  # categorical
    'd_call',  # categorical
    'j_call',  # categorical
    'locus',  # categorical
    'stop',  # bool
    'productive',  # bool
    'UTR',
    'leader',
    'CDR1_nt',
    'CDR1_aa',
    'CDR2_nt',
    'CDR2_aa',
    'cdr3',
    'cdr3_aa',
    'V_nt',
    'V_aa',
    'V_end',
    'np1',
    'D_region',
    'np2',
    'J_nt',
    'sequence_id',
    'barcode',
    'race_G',
    'sequence',
]

_INTEGER_COLUMNS = ('V_errors', 'D_errors', 'J_errors', 'V_CDR3_start')

_RENAME = {
    "CDR3_clusters": "clonotypes",
    "name": "sequence_id",
    "genomic_sequence": "sequence",
    "VD_junction": "np1",
    "DJ_junction": "np2",
    "V_gene": "v_call",
    "D_gene": "d_call",
    "J_gene": "j_call",
    "V_evalue": "v_support",
    "D_evalue": "d_support",
    "J_evalue": "j_support",
    "CDR1_nt": "cdr1",
    "CDR1_aa": "cdr1_aa",
    "CDR2_nt": "cdr2",
    "CDR2_aa": "cdr2_aa",
    "CDR3_nt": "cdr3",
    "CDR3_aa": "cdr3_aa",
}

# The recomputable columns are no longer stored directly in the
# input file, but instead they are recomputed from the other columns
# if requested through usecols
RECOMPUTABLE_COLUMNS = {
    "VDJ_nt": ["v_sequence_start", "j_sequence_end", "sequence"],
    "VDJ_aa": ["v_sequence_start", "j_sequence_end", "sequence"],
    "V_nt": ["v_sequence_alignment"],
    "V_aa": ["v_sequence_alignment_aa"],
    "D_region": ["d_sequence_alignment"],
    "J_nt": ["j_sequence_alignment"],
    "J_aa": ["j_sequence_alignment_aa"],
}


def fix_columns(d, usecols, recompute_cols):
    """
    Changes DataFrame in-place
    """
    # Rename legacy columns
    d.rename(columns=_RENAME, inplace=True)

    # Convert all string columns to str to avoid a PerformanceWarning
    for col in _STRING_COLUMNS:
        if col not in d:
            continue
        d[col].fillna("", inplace=True)
        d[col] = d[col].astype('str')
        # Empty strings have been set to NaN by read_table. Replacing
        # by the empty string avoids problems with groupby, which
        # ignores NaN values.
    # Columns that have any NaN values in them cannot be converted to
    # int due to a numpy limitation.
    for col in _INTEGER_COLUMNS:
        if col not in d.columns:
            continue
        if all(d[col].notnull()):
            d[col] = d[col].astype(int)

    if recompute_cols:
        if "VDJ_nt" in recompute_cols or "VDJ_aa" in recompute_cols:
            vdj_nt = vdj_nt_column(d)
        else:
            vdj_nt = None
        for col in recompute_cols:
            if col == "VDJ_nt":
                d[col] = vdj_nt
            elif col == "VDJ_aa":
                d[col] = vdj_nt.map(nt_to_aa, na_action="ignore")
            elif col == "V_nt":
                d[col] = d["v_sequence_alignment"].str.replace("-", "")
            elif col == "V_aa":
                d[col] = d["v_sequence_alignment_aa"]
            elif col == "D_region":
                d[col] = d["d_sequence_alignment"].str.replace("-", "")
            elif col == "J_nt":
                d[col] = d["j_sequence_alignment"].str.replace("-", "")
            elif col == "J_aa":
                d[col] = d["j_sequence_alignment_aa"]

    if usecols:
        for col in usecols:
            assert col in d.columns, col
        d = d[list(usecols)]  # reorder columns
    return d


def read_table(path, usecols=None, log=False, nrows=None, chunksize=None):
    """
    Read in the table created by the parse subcommand (typically named *.tab)
    """
    if usecols:
        # Adjust requested column names in case the input uses old column names
        available_columns = set(pd.read_table(path, nrows=0).columns)
        new_usecols, recompute_cols = transform_usecols(usecols, available_columns)
    else:
        new_usecols = usecols
        recompute_cols = []

    d = pd.read_table(path, dtype=COLUMN_TYPES, usecols=new_usecols, nrows=nrows)
    d = fix_columns(d, usecols, recompute_cols)
    if log:
        logger.info('%s rows in input table', len(d))
    return d


def read_table_chunks(path, usecols=None, chunksize=1000):
    if usecols:
        # Adjust requested column names in case the input uses old column names
        available_columns = set(pd.read_table(path, nrows=0).columns)
        new_usecols, recompute_cols = transform_usecols(usecols, available_columns)
    else:
        new_usecols = usecols
        recompute_cols = []

    for d in pd.read_table(path, dtype=COLUMN_TYPES, usecols=new_usecols, chunksize=chunksize):
        d = fix_columns(d, usecols, recompute_cols)
        yield d


def transform_usecols(usecols, available_columns):
    recompute_cols = set()
    new_to_old = {v: k for k, v in _RENAME.items()}
    new_usecols = set()
    for col in usecols:
        if col not in available_columns and col in new_to_old:
            logger.info(
                "Requested column %s not found, assuming old file format and using %s instead",
                col, new_to_old[col]
            )
            new_usecols.add(new_to_old[col])
        elif col not in available_columns and col in RECOMPUTABLE_COLUMNS:
            new_usecols.update(RECOMPUTABLE_COLUMNS[col])
            recompute_cols.add(col)
        else:
            new_usecols.add(col)
    return new_usecols, recompute_cols


def vdj_nt_column(table):
    """Return a Series with the nucleotide VDJ sequences"""
    def vdj_nt(row):
        if pd.isna(row.v_call) or pd.isna(row.j_call):
            return np.nan

        return row.sequence[row.v_sequence_start - 1 : row.j_sequence_end]

    return table.apply(vdj_nt, axis=1)
