"""
Find V gene sister sequences shared by multiple libraries.
"""
import logging
from collections import Counter
import pandas as pd

logger = logging.getLogger(__name__)


def add_arguments(parser):
	parser.add_argument('--minimum-frequency', '-n', type=int, metavar='N',
		default=None,
		help='Minimum number of datasets in which sequence must occur (default is no. of files divided by two)')
	parser.add_argument('--minimum-db-diff', '-b', type=int, metavar='DIST', default=1,
		help='Use only sequences that have at least DIST differences to the database sequence. Default: %(default)s')
	parser.add_argument('tables', metavar='DISCOVER.TAB',
		help='Table created by the "discover" command (give at least two)', nargs='+')


def main(args):
	if args.minimum_frequency is None:
		minimum_frequency = max((len(args.tables) + 1) // 2, 2)
	else:
		minimum_frequency = args.minimum_frequency
	logger.info('Minimum frequency set to %s', minimum_frequency)

	# Read in tables
	tables = []
	for path in args.tables:
		table = pd.read_csv(path, sep='\t')
		table = table[table.database_diff >= args.minimum_db_diff]
		table = table.dropna()
		tables.append(table)
		if len(table) == 0:
			logger.warn('Table read from %r is empty after filtering out sequences with database diff >= %s.', path, args.minimum_db_diff)

	# Count V sequence occurrences
	counter = Counter()
	for table in tables:
		counter.update(set(table.consensus))

	# Find most frequent occurrences and print result
	print('count', 'gene', 'database_diff', 'sequence', 'names', sep='\t')
	for sequence, frequency in counter.most_common():
		if frequency < minimum_frequency:
			break
		names = []
		gene = None
		for table in tables:
			matching_rows = table[table.consensus == sequence]
			if matching_rows.empty:
				continue
			names.extend(matching_rows.name)
			if gene is None:
				row = matching_rows.iloc[0]
				gene = row.gene
				database_diff = row.database_diff
				#shm = row['V_SHM']
		print(frequency, gene, database_diff, sequence, *names, sep='\t')
