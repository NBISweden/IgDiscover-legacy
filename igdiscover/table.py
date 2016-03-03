"""
Function for reading the table created by the 'parse' subcommand.
"""
import os.path
import logging
import sqlite3
import pandas as pd
import numpy as np
from tempfile import TemporaryDirectory

logger = logging.getLogger(__name__)

_CACHING = False

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
	'name',
	'barcode',
	'race_G',
	'genomic_sequence',
]

_INTEGER_COLUMNS = ('V_errors', 'J_errors')

def _fix_columns(df):
	"""
	Changes DataFrame in-place
	"""
	# Convert all string columns to str to avoid a PerformanceWarning
	for col in _STRING_COLUMNS:
		df[col].fillna('', inplace=True)
		df[col] = df[col].astype('str')
		# Empty strings have been set to NaN by read_csv. Replacing
		# by the empty string avoids problems with groupby, which
		# ignores NaN values.
	# Columns that have any NaN values in them cannot be converted to
	# int due to a numpy limitation.
	for col in _INTEGER_COLUMNS:
		if all(df[col].notnull()):
				df[col] = df[col].astype(int)


def read_table(path, log=False):
	"""
	Read in the table created by the parse subcommand (typically named *.tab)
	"""
	if _CACHING:
		base, ext = os.path.splitext(path)
		h5path = base + '.h5'

		if not os.path.exists(h5path) or os.path.getmtime(path) > os.path.getmtime(h5path):
			logger.info('(Re-)Creating table cache %s', h5path)
			# Two or more processes may try to re-create the table at the same time.
			# We therefore write the table into a temporary file first and then move
			# it atomically.
			with TemporaryDirectory(dir=os.path.dirname(h5path)) as tempdir:
				temp_h5 = os.path.join(tempdir, 'db.h5')
				df = pd.read_csv(path, sep='\t') #true_values=['yes'], false_values=['no'])
				_fix_columns(df)
				df.to_hdf(temp_h5, 'table', complevel=3, complib='zlib')
				os.rename(temp_h5, h5path)
		d = pd.read_hdf(h5path, 'table')
	else:
		d = pd.read_csv(path, sep='\t')
		_fix_columns(d)
	if log:
		logger.info('%s rows in input table', len(d))
	return d
