"""
Function for reading the table created by the 'parse' subcommand.
"""
import os.path
import logging
import sqlite3
import pandas as pd
from tempfile import TemporaryDirectory

logger = logging.getLogger(__name__)

# Thresholds for filtering
V_GENE_COVERAGE = 90  # at least
J_GENE_COVERAGE = 60  # at least
V_GENE_EVALUE = 1E-3  # at most


def read_table(path, filter=True, log=False):
	"""
	Read in the table created by the parse subcommand. Discard following rows
	if filter is True:
	- no J assigned
	- stop codon found
	- V gene coverage too low
	- J gene coverage too low
	- V gene E-value too high
	"""
	base, ext = os.path.splitext(path)
	h5path = base + '.h5'

	if not os.path.exists(h5path) or os.path.getmtime(path) > os.path.getmtime(h5path):
		logger.info('(Re-)Creating table cache %s', h5path)
		# Two or more processes may try to re-create the table at the same time.
		# We therefore write the table into a temporary file first and then move
		# it atomically.
		with TemporaryDirectory(dir=os.path.dirname(h5path)) as tempdir:
			temp_h5 = os.path.join(tempdir, 'db.h5')
			df = pd.read_csv(path, sep='\t')
			# Convert all mixed-type columns to str to avoid a PerformanceWarning
			for col in df.columns:
				if df[col].dtype.name == 'object':
					df[col] = df[col].astype('str')
			df.to_hdf(temp_h5, 'table')
			os.rename(temp_h5, h5path)

	d = pd.read_hdf(h5path, 'table')
	assert len(d) > 0

	if log: logger.info('%s rows in input table', len(d))

	# Allow old-style %SHM column headers
	d.rename(columns=lambda x: x.replace('%SHM', '_SHM'), inplace=True)
	if not filter:
		return d

	# Both V and J must be assigned
	filtered = d.dropna(subset=('V_gene', 'J_gene'))[:]
	if log: logger.info('%s rows have both V and J assignment', len(filtered))
	filtered['V_gene'] = pd.Categorical(filtered['V_gene'])

	# Filter out sequences that have a stop codon
	filtered = filtered[filtered.stop == 'no']
	if log: logger.info('%s of those do not have a stop codon', len(filtered))

	# Filter out sequences with a too low V gene hit E-value
	filtered = filtered[filtered.V_evalue <= V_GENE_EVALUE]
	if log: logger.info('%s of those have an E-value of at most %s', len(filtered), V_GENE_EVALUE)

	# Filter out sequences with too low V gene coverage
	filtered = filtered[filtered.V_covered >= V_GENE_COVERAGE]
	if log: logger.info('%s of those cover the V gene by at least %s%%', len(filtered), V_GENE_COVERAGE)

	# Filter out sequences with too low J gene coverage
	filtered = filtered[filtered.J_covered >= J_GENE_COVERAGE]
	if log: logger.info('%s of those cover the J gene by at least %s%%', len(filtered), J_GENE_COVERAGE)

	return filtered
