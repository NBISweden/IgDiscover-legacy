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
from .cluster import hamming_single_linkage

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--mismatches', default=1, type=int,
		help='No. of allowed mismatches between CDR3 sequences. Default: %(default)s')
	arg('reftable', help='Reference table with parsed and filtered '
		'IgBLAST results (filtered.tab or filtered.tab)')
	arg('querytable', help='Query table with IgBLAST results (assigned.tab)')


def group_by_cdr3(table, mismatches):
	"""
	Cluster the rows of the table by Hamming distance between
	their CDR3 sequences. Yield (index, group) tuples similar 
	to .groupby().
	"""
	# Cluster all unique CDR3s by Hamming distance
	sequences = list(set(table.CDR3_nt))
	clusters = hamming_single_linkage(sequences, mismatches)

	# Create dict that maps CDR3 sequences to a numeric cluster id
	cluster_ids = dict()
	for cluster_id, cdr3s in enumerate(clusters):
		for cdr3 in cdr3s:
			cluster_ids[cdr3] = cluster_id

	# Assign cluster id to each row
	table['cluster_id'] = table['CDR3_nt'].apply(lambda cdr3: cluster_ids[cdr3])
	for index, group in table.groupby('cluster_id'):
		yield group.drop('cluster_id', axis=1)


def main(args):
	columns = ['count', 'V_gene', 'D_gene', 'J_gene', 'CDR3_nt', 'CDR3_aa',
		'V_errors', 'J_errors', 'barcode', 'VDJ_nt']  # TODO D_errors
	reftable = read_table(args.reftable, usecols=columns)
	reftable = reftable[columns]  # reorder columns
	logger.info('Read reference table with %s rows', len(reftable))
	reftable.insert(6, 'CDR3_length', reftable['CDR3_nt'].apply(len))
	querytable = read_table(args.querytable, usecols=columns + ['name'])
	logger.info('Read query table with %s rows', len(querytable))
	if len(querytable) > len(reftable):
		logger.warning('The reference table is smaller than the '
			'query table! Did you swap query and reference?')

	# Determine set of clonotypes (ignoring CDR3 sequence) to query
	groupby = ('V_gene', 'J_gene', 'CDR3_length')
	query_clonotypes = defaultdict(list)
	for row in querytable.itertuples():
		clonotype = (row.V_gene, row.J_gene, len(row.CDR3_nt))
		query_clonotypes[clonotype].append(row)

	print_header = True
	for clonotype, vjlen_group in reftable.groupby(groupby):
		# (v_gene, j_gene, cdr3_length) = clonotype
		if clonotype not in query_clonotypes:
			continue
		if print_header:
			print(vjlen_group.head(0).to_csv(sep='\t', header=True, index=False))
			print_header = False
		for query_row in query_clonotypes[clonotype]:
			print('# query: {}'.format(query_row.name))
			cdr3 = query_row.CDR3_nt
			is_similar = [
				hamming_distance(cdr3, r.CDR3_nt) <= args.mismatches for r in vjlen_group.itertuples()]
			similar_group = vjlen_group.loc[is_similar, :]
			print(similar_group.to_csv(sep='\t', header=False, index=False))
