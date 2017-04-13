"""
Rudimentary discovery of J genes

Print out the most frequent J sequences in FASTA format.
If one sequence is suffix of another, the longer version is kept.
"""
import logging
from .utils import Merger, merge_overlapping, unique_name
from .table import read_table
from sqt import FastaReader
from sqt.align import edit_distance

logger = logging.getLogger(__name__)


MIN_COUNT = 100


def add_arguments(parser):
	arg = parser.add_argument
	arg('--database', metavar='FASTA',
		help='FASTA file with reference gene sequences')
	arg('table', help='Table with parsed and filtered IgBLAST results')


class SequenceInfo:
	def __init__(self, sequence, count=0, max_count=0, cdr3s=None, v_genes=None):
		self.name = 'name'
		self.sequence = sequence
		self.count = count
		self.max_count = max_count  # an upper bound for the count
		self.cdr3s = cdr3s if cdr3s is not None else set()
		self.v_genes = v_genes if v_genes is not None else set()

	def __repr__(self):
		return 'SequenceInfo({sequence!r}, count={count}, max_count={max_count}, ...)'.format(**vars(self))
	# def __lt__(self, other): ...


class SequenceMerger(Merger):
	"""
	Merge sequences that overlap
	"""
	def __init__(self):
		super().__init__()

	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. If they should not be merged,
		None is returned.
		"""
		m = merge_overlapping(s.sequence, t.sequence)
		if m is not None:
			return SequenceInfo(m, max_count=t.max_count + s.max_count)
		else:
			return None


def main(args):
	if args.database:
		with FastaReader(args.database) as fr:
			database = list(fr)
		logger.info('Read %d sequences from %r', len(database), args.database)
	else:
		database = None
	table = read_table(args.table,
		usecols=['count', 'V_gene', 'J_gene', 'V_errors', 'J_nt', 'D_region', 'CDR3_nt', 'genomic_sequence'])
	logger.info('Table with %s rows read', len(table))
	table = table[table.V_errors == 0]
	logger.info('Keeping %s rows that have zero V mismatches', len(table))

	# Merge candidate sequences that overlap. If one candidate is longer than
	# another, this is typically a sign that IgBLAST has not extended the
	# alignment long enough.
	merger = SequenceMerger()
	for sequence, group in table.groupby('J_nt'):
		merger.add(SequenceInfo(sequence, max_count=len(group)))
	logger.info('After merging overlapping sequences, %s remain', len(merger))

	# Use only records that have a chance of reaching the required MIN_COUNT
	records = {info.sequence: info for info in merger if info.max_count >= MIN_COUNT}

	# Speed up search by looking for most common sequences first
	search_order = sorted(records, key=lambda s: records[s].max_count, reverse=True)

	# Count full-text occurrences in the genomic_sequence, circumventing
	# inaccurate IgBLAST alignment boundaries
	for row in table.itertuples():
		# print('current row', row)
		for needle in search_order:
			# print('current needle', needle)
			if needle in row.genomic_sequence:
				record = records[needle]
				record.count += 1
				record.v_genes.add(row.V_gene)
				record.cdr3s.add(row.CDR3_nt)
				break

	# Print output table
	print('name', 'count', 'V_genes', 'CDR3s', 'database', 'database_diff', 'sequence', sep='\t')
	i = 0
	for record in sorted(records.values(), key=lambda r: (r.count, r.sequence), reverse=True):
		# TODO
		# if record.count < MIN_COUNT:
		# 	break
		if database:
			distances = [(edit_distance(d.sequence, record.sequence), d) for d in database]
			distance, closest = min(distances, key=lambda x: x[0])
			db_name = closest.name
			if distance == 0:
				name = closest.name
			else:
				name = unique_name(closest.name, record.sequence)
		else:
			db_name = ''
			distance = -1
			name = unique_name('J', record.sequence)
		print(name,
			record.count,
			len(set(record.v_genes)),
			len(set(record.cdr3s)),
			db_name,
			distance,
			record.sequence,
			sep='\t')
		i += 1
