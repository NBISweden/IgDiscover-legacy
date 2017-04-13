"""
Discover D and J genes

The most frequent D/J sequences are considered candidates.

For J genes, candidate sequences are merged if they overlap each other.
"""
import logging
import pandas as pd
from .utils import Merger, merge_overlapping, unique_name
from .table import read_table
from sqt import FastaReader
from sqt.align import edit_distance

logger = logging.getLogger(__name__)


MIN_COUNT = 100
MINIMUM_CANDIDATE_LENGTH = 5


def add_arguments(parser):
	arg = parser.add_argument
	arg('--database', metavar='FASTA',
		help='FASTA file with reference gene sequences')
	arg('--merge', default=None, action='store_true', help='Merge overlapping genes. '
		'Default: Enabled for J, disabled for D.')
	arg('--no-merge', dest='merge', action='store_false', help='Do not merge overlapping genes')
	arg('--gene', default='J', choices=('D', 'J'),
		help='Which gene category to discover. Default: %(default)s')
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
		usecols=['count', 'V_gene', 'J_gene', 'V_errors', 'J_nt', 'D_region', 'CDR3_nt'])
	logger.info('Table with %s rows read', len(table))
	table = table[table.V_errors == 0]
	logger.info('Keeping %s rows that have zero V mismatches', len(table))

	if args.merge is None:
		args.merge = args.gene == 'J'

	column = 'J_nt' if args.gene == 'J' else 'D_region'
	if args.merge:
		# Merge candidate sequences that overlap. If one candidate is longer than
		# another, this is typically a sign that IgBLAST has not extended the
		# alignment long enough.
		merger = SequenceMerger()
		for sequence, group in table.groupby(column):
			if len(sequence) < MINIMUM_CANDIDATE_LENGTH:
				continue
			merger.add(SequenceInfo(sequence, max_count=len(group)))
		logger.info('After merging overlapping %s sequences, %s remain', args.gene, len(merger))
		candidates = list(merger)
	else:
		candidates = []
		for sequence, group in table.groupby(column):
			if len(sequence) < MINIMUM_CANDIDATE_LENGTH:
				continue
			candidates.append(SequenceInfo(sequence, max_count=len(group)))
		logger.info('Collected %s unique %s sequences', len(candidates), args.gene)

	# Use only records that have a chance of reaching the required MIN_COUNT
	records = {info.sequence: info for info in candidates if info.max_count >= MIN_COUNT}

	# Speed up search by looking for most common sequences first
	search_order = sorted(records, key=lambda s: records[s].max_count, reverse=True)

	# Count full-text occurrences in the genomic_sequence, circumventing
	# inaccurate IgBLAST alignment boundaries
	# TODO limit the search to the gene region (especially for D genes)
	del table
	cols = ['V_gene', 'V_errors', 'CDR3_nt', 'genomic_sequence']
	for chunk in pd.read_csv(args.table, usecols=cols, chunksize=10000, sep='\t'):
		chunk = chunk[chunk.V_errors == 0]
		for row in chunk.itertuples():
			for needle in search_order:
				if needle in row.genomic_sequence:
					record = records[needle]
					record.count += 1
					record.v_genes.add(row.V_gene)
					record.cdr3s.add(row.CDR3_nt)
					if args.merge:
						break

	# Print output table
	print('name', 'count', 'V_genes', 'CDR3s', 'database', 'database_diff', 'sequence', sep='\t')
	i = 0
	for record in sorted(records.values(), key=lambda r: (r.count, r.sequence), reverse=True):
		if record.count < 1:
			break
		if database:
			distances = [(edit_distance(db.sequence, record.sequence), db) for db in database]
			distance, closest = min(distances, key=lambda x: x[0])
			db_name = closest.name
			if distance == 0:
				name = closest.name
			else:
				name = unique_name(closest.name, record.sequence)
		else:
			db_name = ''
			distance = -1
			name = unique_name(args.gene, record.sequence)
		print(name,
			record.count,
			len(set(record.v_genes)),
			len(set(record.cdr3s)),
			db_name,
			distance,
			record.sequence,
			sep='\t')
		i += 1
