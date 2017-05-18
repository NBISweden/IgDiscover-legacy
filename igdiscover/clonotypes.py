"""
Group a table of assigned sequences by clonotype

Two sequences have the same clonotype if
- their V and J assignments are the same
- the length of their CDR3 is identical
- the difference between their CDR3s (in terms of mismatches)
  is not higher than a given threshold (by default 1)

The output is a table with one row per input sequence, sorted by
clonotype, plus an empty line between each group of sequences
that have the same clonotype.

The table is written to standard output.
"""
import logging

from .table import read_table
from .cluster import hamming_single_linkage

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--mismatches', default=1, type=int,
		help='No. of allowed mismatches between CDR3 sequences. Default: %(default)s')
	arg('table', help='Table with parsed and filtered IgBLAST results')


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
	table = read_table(args.table, usecols=columns)
	logger.info('Read table with %s rows', len(table))
	table.insert(6, 'CDR3_length', table['CDR3_nt'].apply(len))
	print_header = True
	prev_v = None
	n = 0
	for (v_gene, j_gene, cdr3_length), vj_group in table.groupby(('V_gene', 'J_gene', 'CDR3_length')):
		if prev_v != v_gene:
			logger.info('%d clonotypes seen. Processing %s', n, v_gene)
		prev_v = v_gene
		for group in group_by_cdr3(vj_group.copy(), mismatches=args.mismatches):
			# We get an intentional empty line between groups since
			# to_csv() already includes a line break
			print(group.to_csv(sep='\t', header=print_header, index=False))
			print_header = False
			n += 1
	logger.info('A total of %d clonotypes was seen', n)
