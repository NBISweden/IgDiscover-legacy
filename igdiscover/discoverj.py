"""
Discover D and J genes

The most frequent D/J sequences are considered candidates.

For J genes, candidate sequences are merged if they overlap each other.

The result table is written to standard output. Use --fasta to also
generate FASTA output.
"""
import logging
import pandas as pd
from sqt import FastaReader
from sqt.align import edit_distance
from .utils import Merger, merge_overlapping, unique_name, is_same_gene
from .table import read_table

logger = logging.getLogger(__name__)

MIN_COUNT = 100
MINIMUM_CANDIDATE_LENGTH = 5


def add_arguments(parser):
	arg = parser.add_argument
	arg('--database', metavar='FASTA',
		help='FASTA file with reference gene sequences')
	arg('--merge', default=None, action='store_true', help='Merge overlapping genes. '
		'Default: Enabled for J, disabled for D and V.')
	arg('--no-merge', dest='merge', action='store_false', help='Do not merge overlapping genes')
	arg('--gene', default='J', choices=('V', 'D', 'J'),
		help='Which gene category to discover. Default: %(default)s')
	arg('--allele-ratio', type=float, metavar='RATIO', default=0.1,
		help='Required allele ratio. Works only for genes named "NAME*ALLELE". Default: %(default)s')
	arg('--cross-mapping-ratio', type=float, metavar='RATIO', default=None,
		help='Ratio for detection of cross-mapping artifacts. Default: %(default)s')
	arg('--fasta', help='Write discovered sequences to FASTA file')
	arg('table', help='Table with parsed and filtered IgBLAST results')


class SequenceInfo:
	__slots__ = ('name', 'sequence', 'count', 'max_count', 'cdr3s', 'other_genes', 'db_name',
		'db_distance')

	def __init__(self, name, sequence, count=0, max_count=0, cdr3s=None, other_genes=None,
			db_name=None, db_distance=None):
		self.name = name
		self.sequence = sequence
		self.count = count
		self.max_count = max_count  # an upper bound for the count
		self.cdr3s = cdr3s if cdr3s is not None else set()
		self.other_genes = other_genes if other_genes is not None else set()
		self.db_name = db_name
		self.db_distance = db_distance

	def __repr__(self):
		return 'SequenceInfo({sequence!r}, count={count}, max_count={max_count}, ...)'.format(**vars(self))
	# def __lt__(self, other): ...


class OverlappingSequenceMerger(Merger):
	"""
	Merge sequences that overlap
	"""
	def __init__(self, cross_mapping_ratio=None):
		super().__init__()
		self._cross_mapping_ratio = cross_mapping_ratio

	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. If they should not be merged,
		None is returned.
		"""
		m = merge_overlapping(s.sequence, t.sequence)
		if m is not None:
			return SequenceInfo(s.name, m, max_count=t.max_count + s.max_count)

		return None


class AlleleRatioMerger(Merger):
	"""
	Discard sequences with too low allele ratio
	"""
	def __init__(self, allele_ratio, cross_mapping_ratio):
		super().__init__()
		self._allele_ratio = allele_ratio
		self._cross_mapping_ratio = cross_mapping_ratio

	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. If they should not be merged,
		None is returned.
		"""
		# TODO copy-and-pasted from germlinefilter
		#
		# Check allele ratio. Somewhat similar to cross-mapping, but
		# this uses sequence names to decide whether two genes can be
		# alleles of each other and the ratio is between the CDR3s_exact
		# values
		if self._allele_ratio and is_same_gene(s.name, t.name):
			for u, v in [(s, t), (t, s)]:
				ratio = len(set(u.cdr3s)) / len(set(v.cdr3s))
				if ratio < self._allele_ratio:
					# logger.info('Allele ratio %.4f too low for %r compared to %r',
					# 	ratio, u.name, v.name)
					return v

		if self._cross_mapping_ratio:
			# When checking for cross mapping, ignore overhanging bases in the 5' end.
			# Example:
			# ---ACTACGACTA...
			# XXX|||||X||||
			# ATTACTACTACTA...
			if len(t.sequence) < len(s.sequence):
				t, s = s, t  # s is now the shorter sequence
			t_seq = t.sequence[len(t.sequence) - len(s.sequence):]
			s_seq = s.sequence
			dist = edit_distance(s_seq, t_seq, 1)
			if dist > 1:
				return None
			total_count = (s.count + t.count)
			if total_count == 0:
				return None
			for u, v in [(s, t), (t, s)]:
				ratio = u.count / total_count
				if ratio < self._cross_mapping_ratio:
					# u is probably a cross-mapping artifact of the higher-expressed v
					logger.info('%r is a cross-mapping artifact of %r (ratio %.4f)',
						u.name, v.name, ratio)
					return v

		return None


def count_full_text_occurrences(candidates, table_path, other_gene, other_errors, merge):
	# Use only records that have a chance of reaching the required MIN_COUNT
	records = {info.sequence: info for info in candidates if info.max_count >= MIN_COUNT}

	# Count full-text occurrences in the genomic_sequence, circumventing
	# inaccurate IgBLAST alignment boundaries
	# TODO limit the search to the gene region (especially for D genes)
	# Speed up search by looking for most common sequences first
	search_order = sorted(records, key=lambda s: records[s].max_count, reverse=True)
	cols = [other_gene, 'V_errors', 'J_errors', 'CDR3_nt', 'genomic_sequence']
	for chunk in pd.read_csv(table_path, usecols=cols, chunksize=10000, sep='\t'):
		chunk = chunk[chunk[other_errors] == 0]
		for row in chunk.itertuples():
			for needle in search_order:
				if needle in row.genomic_sequence:
					record = records[needle]
					record.count += 1
					record.other_genes.add(getattr(row, other_gene))
					record.cdr3s.add(row.CDR3_nt)
					if merge:
						break
	return records.values()


def print_table(records, other_gene):
	# Print output table
	print('name', 'count', other_gene + 's', 'CDR3s', 'database', 'database_diff', 'sequence', sep='\t')
	for record in records:
		print(record.name,
			record.count,
			len(set(record.other_genes)),
			len(set(record.cdr3s)),
			record.db_name if record.db_name is not None else '',
			record.db_distance if record.db_distance is not None else -1,
			record.sequence,
			sep='\t')


def main(args):
	if args.database:
		with FastaReader(args.database) as fr:
			database = list(fr)
		logger.info('Read %d sequences from %r', len(database), args.database)
	else:
		database = None
	if args.gene == 'V':
		column = 'V_nt'
	elif args.gene == 'J':
		column = 'J_nt'
	else:
		assert args.gene == 'D'
		column = 'D_region'
	other = 'V' if args.gene in ('D', 'J') else 'J'
	other_gene = other + '_gene'
	other_errors = other + '_errors'
	table = read_table(args.table,
		usecols=['count', 'V_gene', 'J_gene', 'V_errors', 'J_errors', column, 'CDR3_nt'])
	logger.info('Table with %s rows read', len(table))
	table = table[table[other_errors] == 0]
	logger.info('Keeping %s rows that have no %s mismatches', len(table), other)

	if args.merge is None:
		args.merge = args.gene == 'J'

	if args.merge:
		logger.info('Merging overlapping sequences ...')
		# Merge candidate sequences that overlap. If one candidate is longer than
		# another, this is typically a sign that IgBLAST has not extended the
		# alignment long enough.
		merger = OverlappingSequenceMerger()
		for sequence, group in table.groupby(column):
			if len(sequence) < MINIMUM_CANDIDATE_LENGTH:
				continue
			merger.add(SequenceInfo(None, sequence, max_count=len(group)))
		logger.info('After merging overlapping %s sequences, %s remain', args.gene, len(merger))
		candidates = list(merger)
	else:
		candidates = []
		for sequence, group in table.groupby(column):
			if len(sequence) < MINIMUM_CANDIDATE_LENGTH:
				continue
			candidates.append(SequenceInfo(None, sequence, max_count=len(group)))
		logger.info('Collected %s unique %s sequences', len(candidates), args.gene)
	del table

	logger.info('Counting occurrences ...')
	records = count_full_text_occurrences(candidates, args.table, other_gene, other_errors, args.merge)

	# Assign names etc.
	if database:
		for record in records:
			distances = [(edit_distance(db.sequence, record.sequence), db) for db in database]
			record.db_distance, closest = min(distances, key=lambda x: x[0])
			record.db_name = closest.name
			if record.db_distance == 0:
				record.name = closest.name
			else:
				record.name = unique_name(closest.name, record.sequence)
		# Merge by allele ratio
	else:
		for record in records:
			record.name = unique_name(args.gene, record.sequence)

	if args.allele_ratio or args.cross_mapping_ratio:
		arm = AlleleRatioMerger(args.allele_ratio, args.cross_mapping_ratio)
		arm.extend(records)
		records = list(arm)
		logger.info('After filtering by allele ratio and/or cross-mapping ratio, %d candidates remain', len(records))

	records = sorted(records, key=lambda r: (r.count, r.sequence), reverse=True)
	records = [r for r in records if r.count > 1]

	print_table(records, other_gene)
	if args.fasta:
		with open(args.fasta, 'w') as f:
			for record in sorted(records, key=lambda r: r.name):
				print('>{}\n{}'.format(record.name, record.sequence), file=f)

	logger.info('Wrote %d genes', len(records))
