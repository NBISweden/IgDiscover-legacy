"""
Query two ...
"""
import logging
from collections import defaultdict
from contextlib import ExitStack

import pandas as pd
from xopen import xopen

from igdiscover.table import read_table
from igdiscover.utils import slice_arg
from igdiscover.clonotypes import is_similar_with_junction, CLONOTYPE_COLUMNS

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	# arg('--minimum-count', '-c', metavar='N', default=1, type=int,
	# 	help='Discard all rows with count less than N. Default: %(default)s')
	# arg('--cdr3-core', default=None,
	# 	type=slice_arg, metavar='START:END',
	# 	help='START:END defines the non-junction region of CDR3 '
	# 		'sequences. Use negative numbers for END to count '
	# 		'from the end. Regions before and after are considered to '
	# 		'be junction sequence, and for two CDR3s to be considered '
	# 		'similar, at least one of the junctions must be identical. '
	# 		'Default: no junction region.')
	# arg('--mismatches', default=1, type=float,
	# 	help='No. of allowed mismatches between CDR3 sequences. '
	# 		'Can also be a fraction between 0 and 1 (such as 0.15), '
	# 		'interpreted relative to the length of the CDR3 (minus the front non-core). '
	# 		'Default: %(default)s')
	arg('--aa', default=False, action='store_true',
		help='Count CDR3 mismatches on amino-acid level. Default: Compare nucleotides.')
	arg('tables', nargs='+', help='clonotype tables')


def collect(querytable, reftable, mismatches, cdr3_core_slice, cdr3_column):
	"""
	Find all queries from the querytable in the reftable.

	Yield tuples (query_rows, similar_rows) where the query_rows is a list
	with all the rows that have the same result. similar_rows is a DataFrame
	whose rows are the ones matching the query.
	"""
	# The vjlentype is a "clonotype without CDR3 sequence" (only V, J, CDR3 length)
	# Determine set of vjlentypes to query
	query_vjlentypes = defaultdict(list)
	for row in querytable.itertuples():
		vjlentype = (row.V_gene, row.J_gene, len(row.CDR3_nt))
		query_vjlentypes[vjlentype].append(row)

	groupby = ('V_gene', 'J_gene', 'CDR3_length')
	for vjlentype, vjlen_group in reftable.groupby(groupby):
		# (v_gene, j_gene, cdr3_length) = vjlentype
		if vjlentype not in query_vjlentypes:
			continue

		# Collect results for this vjlentype. The results dict
		# maps row indices (into the vjlen_group) to each query_row,
		# allowing us to group identical results together.
		results = defaultdict(list)
		for query_row in query_vjlentypes.pop(vjlentype):
			cdr3 = getattr(query_row, cdr3_column)
			# Save indices of the rows that are similar to this query
			indices = tuple(index for index, r in enumerate(vjlen_group.itertuples())
				if is_similar_with_junction(cdr3, getattr(r, cdr3_column), mismatches, cdr3_core_slice))
			results[indices].append(query_row)

		# Yield results, grouping queries that lead to the same result
		for indices, query_rows in results.items():
			if not indices:
				for query_row in query_rows:
					yield ([query_row], [])
				continue

			similar_group = vjlen_group.iloc[list(indices), :]
			yield (query_rows, similar_group)

	# Yield result tuples for all the queries that have not been found
	for queries in query_vjlentypes.values():
		for query_row in queries:
			yield ([query_row], [])


def main(args):
	cdr3_column = 'CDR3_aa' if args.aa else 'CDR3_nt'
	# CDR3_length is unused, but this ensures that we are reading the right type of file
	usecols = ['count', 'V_gene', 'J_gene', 'CDR3_nt', 'CDR3_length', cdr3_column]
	tables = []
	for path in args.tables:
		table = read_table(path, usecols=usecols).set_index(
			['V_gene', 'J_gene', cdr3_column])
		tables.append(table)
		logger.info('Read table from %r with %s rows', path, len(table))

	common = tables[0].join(tables[1], lsuffix='_left', rsuffix='_right', how='inner')
	print('Common clonotypes:', len(common))

	n = common[['count_left', 'count_right']].min(axis=1).sum()
	print('Common sequences:', n)


if __name__ == '__main__':
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
