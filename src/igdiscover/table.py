"""
Function for reading the table created by the 'parse' subcommand.
"""
import logging
import pandas as pd

logger = logging.getLogger(__name__)


# These columns contain string data
# convert them to str to avoid a PerformanceWarning
# TODO some of these are actually categorical or bool
_STRING_COLUMNS = [
    'V_gene',  # categorical
    'D_gene',  # categorical
    'J_gene',  # categorical
    'chain',  # categorical
    'stop',  # bool
    'productive',  # bool
    'UTR',
    'leader',
    'CDR1_nt',
    'CDR1_aa',
    'CDR2_nt',
    'CDR2_aa',
    'CDR3_nt',
    'CDR3_aa',
    'V_nt',
    'V_aa',
    'V_end',
    'VD_junction',
    'D_region',
    'DJ_junction',
    'J_nt',
    'VDJ_nt',
    'VDJ_aa',
    'name',
    'barcode',
    'race_G',
    'genomic_sequence',
]

_INTEGER_COLUMNS = ('V_errors', 'D_errors', 'J_errors', 'V_CDR3_start')


def fix_columns(df):
    """
    Changes DataFrame in-place
    """
    # Convert all string columns to str to avoid a PerformanceWarning
    for col in _STRING_COLUMNS:
        if col not in df:
            continue
        df[col].fillna('', inplace=True)
        df[col] = df[col].astype('str')
        # Empty strings have been set to NaN by read_csv. Replacing
        # by the empty string avoids problems with groupby, which
        # ignores NaN values.
    # Columns that have any NaN values in them cannot be converted to
    # int due to a numpy limitation.
    for col in _INTEGER_COLUMNS:
        if col not in df.columns:
            continue
        if all(df[col].notnull()):
            df[col] = df[col].astype(int)

    # TODO backwards compatibility
    if 'CDR3_clusters' in df.columns:
        df.rename(columns={'CDR3_clusters': 'clonotypes'}, inplace=True)


def read_table(path, usecols=None, log=False, nrows=None):
    """
    Read in the table created by the parse subcommand (typically named *.tab)
    """
    d = pd.read_csv(path, usecols=usecols, sep='\t', nrows=nrows)
    fix_columns(d)
    if log:
        logger.info('%s rows in input table', len(d))
    return d
