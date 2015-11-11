"""
Find common V genes between two different antibody libraries.
"""
import logging
from collections import Counter

from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	parser.add_argument('--minimum-frequency', '-n', type=int, metavar='N',
		default=None,
		help='Minimum number of datasets in which sequence must occur (default is no. of files divided by two)')
	parser.add_argument('table', help='Tables with parsed and filtered IgBLAST results (give at least two)', nargs='+')


def main(args):
	if args.minimum_frequency is None:
		# args.table is a list of file names
		minimum_frequency = max((len(args.table) + 1) // 2, 2)
	else:
		minimum_frequency = args.minimum_frequency
	logger.info('Minimum frequency set to %s', minimum_frequency)

	# Read in tables
	tables = []
	for path in args.table:
		table = read_table(path)
		table = table.loc[:,['V_gene', 'V_SHM', 'V_nt', 'name']]
		tables.append(table)

	# Count V sequence occurrences
	counter = Counter()
	for table in tables:
		counter.update(set(table.V_nt))

	# Find most frequent occurrences and print result
	print('Frequency', 'Gene', '%SHM', 'Sequence', sep='\t')
	for sequence, frequency in counter.most_common():
		if frequency < minimum_frequency:
			break
		names = []
		gene = None
		for table in tables:
			matching_rows = table[table.V_nt == sequence]
			if matching_rows.empty:
				continue
			names.extend(matching_rows.name)
			if gene is None:
				row = matching_rows.iloc[0]
				gene = row['V_gene']
				shm = row['V_SHM']
		print(frequency, gene, shm, sequence, *names, sep='\t')
