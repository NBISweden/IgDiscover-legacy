"""
Query a table of assigned sequences by clonotype

Two sequences have the same clonotype if
- their V and J assignments are the same
- the length of their CDR3 is identical
- the difference between their CDR3s (in terms of mismatches)
  is not higher than a given threshold (by default 1)

Clonotypes for the query sequences are determined and sequences
in the input table that have this clonotype are reported.

The table is written to standard output.
"""
import logging
from collections import defaultdict
from sqt.align import hamming_distance

from .table import read_table
from .utils import slice_arg

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--minimum-count', '-c', metavar='N', default=1, type=int,
		help='Discard all rows with count less than N. Default: %(default)s')
	arg('--cdr3-core', default=None,
		type=slice_arg, metavar='START:END',
		help='START:END defines the non-junction region of CDR3 '
			'sequences. Use negative numbers for END to count '
			'from the end. Regions before and after are considered to '
			'be junction sequence, and for two CDR3s to be considered '
			'similar, at least one of the junctions must be identical. '
			'Default: no junction region.')
	arg('--mismatches', default=1, type=int,
		help='No. of allowed mismatches between CDR3 sequences. Default: %(default)s')
	arg('reftable', help='Reference table with parsed and filtered '
		'IgBLAST results (filtered.tab)')
	arg('querytable', help='Query table with IgBLAST results (assigned.tab or filtered.tab)')


def is_similar_with_junction(s, t, mismatches, cdr3_core):
	"""
	Return whether strings s and t have at most the given number of mismatches
	*and* have at least one identical junction.
	"""
	distance_ok = hamming_distance(s, t) <= mismatches
	if cdr3_core is None:
		return distance_ok
	return distance_ok and (
			(s[:cdr3_core.start] == t[:cdr3_core.start]) or
			(s[cdr3_core.stop:] == t[cdr3_core.stop:]))


def main(args):
	columns = ['name', 'count', 'V_gene', 'D_gene', 'J_gene', 'CDR3_nt', 'CDR3_aa',
		'V_errors', 'J_errors', 'V_SHM', 'J_SHM', 'barcode', 'VDJ_nt', 'VDJ_aa']  # TODO D_errors
	querytable = read_table(args.querytable, usecols=columns)
	querytable = querytable[columns]  # reorder columns
	logger.info('Read query table with %s rows', len(querytable))
	reftable = read_table(args.reftable, usecols=columns)
	reftable = reftable[columns]
	logger.info('Read reference table with %s rows', len(reftable))
	if args.minimum_count > 1:
		reftable = reftable[reftable['count'] >= args.minimum_count]
		logger.info('After filtering out rows with count < %s, %s rows remain', args.minimum_count,
			len(reftable))
	for tab in querytable, reftable:
		tab.insert(6, 'CDR3_length', tab['CDR3_nt'].apply(len))

	if len(querytable) > len(reftable):
		logger.warning('The reference table is smaller than the '
			'query table! Did you swap query and reference?')

	# The vjlentype is a "clonotype without CDR3 sequence" (only V, J, CDR3 length)
	# Determine set of vjlentypes to query
	groupby = ('V_gene', 'J_gene', 'CDR3_length')
	query_vjlentypes = defaultdict(list)
	for row in querytable.itertuples():
		vjlentype = (row.V_gene, row.J_gene, len(row.CDR3_nt))
		query_vjlentypes[vjlentype].append(row)

	print_header = True
	for vjlentype, vjlen_group in reftable.groupby(groupby):
		# (v_gene, j_gene, cdr3_length) = vjlentype
		if vjlentype not in query_vjlentypes:
			continue
		if print_header:
			if args.cdr3_core:
				print('keep?', end='\t')
			print(vjlen_group.head(0).to_csv(sep='\t', header=True, index=False), sep='\t')
			print_header = False

		# Collect results for this vjlentype
		results = defaultdict(list)
		for query_row in query_vjlentypes.pop(vjlentype):
			cdr3 = query_row.CDR3_nt
			# TODO use is_similar_with_junction
			# Save indices of the rows that are similar to this query
			indices = tuple(index for index, r in enumerate(vjlen_group.itertuples())
				if hamming_distance(cdr3, r.CDR3_nt) <= args.mismatches)
			results[indices].append(query_row)

		# Print results, grouping queries that lead to the same result
		for indices, query_rows in results.items():
			for query_row in query_rows:
				if args.cdr3_core:
					print('', end='\t')
				print('# Query: {}'.format(query_row.name), '', *(query_row[3:]), sep='\t')
				if not indices:
					print()

			similar_group = vjlen_group.iloc[indices, :]
			if args.cdr3_core is not None:
				cdr3_core = args.cdr3_core
				cdr3_head = cdr3[:cdr3_core.start]
				cdr3_tail = cdr3[cdr3_core.stop:]
				similar_group = similar_group.copy()
				similar_group.insert(0, 'junction_ident', [
					int((cdr3_head == r.CDR3_nt[:cdr3_core.start]) or
					(cdr3_tail == r.CDR3_nt[cdr3_core.stop:])) for r in similar_group.itertuples()])
			if indices:
				print(similar_group.to_csv(sep='\t', header=False, index=False))

	# Output all the queries that have not been found
	for queries in query_vjlentypes.values():
		for query_row in queries:
			if args.cdr3_core:
				print('', end='\t')
			print('# Query: {}'.format(query_row.name), '', *(query_row[3:]), sep='\t', end='\n\n')
