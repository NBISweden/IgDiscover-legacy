"""
Function for reading the table created by the 'parse' subcommand.
"""
import os.path
import logging
import sqlite3
import pandas as pd

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
	sdb = base + '.sdb'
	if not os.path.exists(sdb) or os.path.getmtime(path) > os.path.getmtime(sdb):
		# SQLite3 file missing or out of date
		logger.info('(Re-)Creating table cache %s', sdb)
		df = pd.read_csv(path, sep='\t')
		# Delete the file if it exists. This is faster than using
		# if_exists='replace' in the to_sql call below.
		try:
			os.remove(sdb)
		except FileNotFoundError:
			pass
		connection = sqlite3.connect(sdb)
		df.to_sql('data', connection)
		del df
		connection.commit()
	else:
		connection = sqlite3.connect(sdb)

	d = pd.read_sql('SELECT * FROM data', connection)

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
