"""
Function for reading the table created by the 'parse' subcommand.
"""
import logging

import numpy as np
import pandas as pd

from igdiscover.utils import nt_to_aa

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
    'chain',  # categorical
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


def fix_columns(df):
    """
    Changes DataFrame in-place
    """
    # Rename legacy columns
    df.rename(columns=_RENAME, inplace=True)

    # Convert all string columns to str to avoid a PerformanceWarning
    for col in _STRING_COLUMNS:
        if col not in df:
            continue
        df[col].fillna('', inplace=True)
        df[col] = df[col].astype('str')
        # Empty strings have been set to NaN by read_table. Replacing
        # by the empty string avoids problems with groupby, which
        # ignores NaN values.
    # Columns that have any NaN values in them cannot be converted to
    # int due to a numpy limitation.
    for col in _INTEGER_COLUMNS:
        if col not in df.columns:
            continue
        if all(df[col].notnull()):
            df[col] = df[col].astype(int)


def read_table(path, usecols=None, log=False, nrows=None):
    """
    Read in the table created by the parse subcommand (typically named *.tab)
    """
    recompute_cols = []
    new_usecols = usecols
    if usecols:
        # Adjust requested column names in case the input uses old column names
        available_columns = set(pd.read_table(path, sep="\t", nrows=0).columns)

        new_to_old = {v: k for k, v in _RENAME.items()}
        new_usecols = []
        for col in usecols:
            if col not in available_columns and col in new_to_old:
                logger.info(
                    "Requested column %s not found, assuming old file format and using %s instead",
                    col, new_to_old[col]
                )
                new_usecols.append(new_to_old[col])
            else:
                new_usecols.append(col)

        if "VDJ_nt" in usecols or "VDJ_aa" in usecols:
            # We no longer store these columns in the file, derive from other columns
            for col in "VDJ_nt", "VDJ_aa":
                if col in new_usecols:
                    new_usecols.remove(col)
                    recompute_cols.append(col)
            new_usecols.extend(["v_sequence_start", "j_sequence_end", "sequence"])

    d = pd.read_table(path, dtype=COLUMN_TYPES, usecols=new_usecols, nrows=nrows)
    fix_columns(d)

    if recompute_cols:
        vdj_nt = vdj_nt_column(d)
        for col in recompute_cols:
            if col == "VDJ_nt":
                d[col] = vdj_nt
            elif col == "VDJ_aa":
                d[col] = vdj_nt.map(nt_to_aa, na_action="ignore")

    if usecols:
        for col in usecols:
            assert col in d.columns, col
        d = d[list(usecols)]  # reorder columns
    if log:
        logger.info('%s rows in input table', len(d))
    return d


def vdj_nt_column(table):
    """Return a Series with the nucleotide VDJ sequences"""
    def vdj_nt(row):
        if pd.isna(row.v_call) or pd.isna(row.j_call):
            return np.nan

        return row.sequence[row.v_sequence_start - 1 : row.j_sequence_end]

    return table.apply(vdj_nt, axis=1)
