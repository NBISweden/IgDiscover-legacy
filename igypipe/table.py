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
	tab_timestamp = os.path.getmtime(path)
	base, ext = os.path.splitext(path)
	sdb = base + '.sdb'
	# Long timeout since another process may just be creating the table
	logger.info('Trying to connect...')
	connection = sqlite3.connect(sdb, isolation_level='exclusive', timeout=120)
	connection.set_trace_callback(logger.info)
	logger.info('Connected successfully')
	connection.execute('BEGIN EXCLUSIVE TRANSACTION')
	logger.info('started exclusive transaction')

	try:
		rows = list(connection.execute('SELECT stamp FROM timestamps'))
		logger.info('successfully read from timestamps')
		out_of_date = len(rows) < 1 or rows[0][0] < tab_timestamp
		logger.info('timestamp in tab: %s. of file: %s. outofdate: %s', rows[0][0], tab_timestamp, out_of_date)
	except sqlite3.OperationalError:
		# Table missing
		logger.info('table missing')
		out_of_date = True
	if out_of_date:
		logger.info('(Re-)Creating table cache %s', sdb)
		import time; time.sleep(10)
		logger.info('sleeping done')
		connection.execute('CREATE TABLE IF NOT EXISTS timestamps (stamp FLOAT)')
		connection.execute('DELETE FROM timestamps')
		connection.execute('INSERT INTO timestamps VALUES (?)', (tab_timestamp,))
		logger.info('pd.read_csv')
		df = pd.read_csv(path, sep='\t')
		logger.info('df.to_sql')
		df.to_sql('data', connection, if_exists='replace')
		logger.info('connection.commit')
		connection.execute('COMMIT')#commit()
		del df
	else:
		logger.info('table cache up to date')
	d = pd.read_sql('SELECT * FROM data', connection)
	connection.close()

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
