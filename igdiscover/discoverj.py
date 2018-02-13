"""
Discover D and J genes

The most frequent D/J sequences are considered candidates.

For J genes, candidate sequences are merged if they overlap each other.

The result table is written to standard output. Use --fasta to also
generate FASTA output.
"""
import logging
import functools
from collections import Counter
import pandas as pd
from sqt import FastaReader
from sqt.align import edit_distance
from .utils import Merger, merge_overlapping, unique_name, is_same_gene, slice_arg
from .table import read_table, fix_columns

logger = logging.getLogger(__name__)

MINIMUM_CANDIDATE_LENGTH = 5


def add_arguments(parser):
	arg = parser.add_argument
	arg('--database', metavar='FASTA',
		help='FASTA file with reference gene sequences')
	arg('--merge', default=None, action='store_true', help='Merge overlapping genes. '
		'Default: Enabled for D, disabled for J and V.')
	arg('--no-merge', dest='merge', action='store_false', help='Do not merge overlapping genes')
	arg('--gene', default='J', choices=('V', 'D', 'J'),
		help='Which gene category to discover. Default: %(default)s')
	arg('--allele-ratio', type=float, metavar='RATIO', default=0.1,
		help='Required allele ratio. Works only for genes named "NAME*ALLELE". Default: %(default)s')
	arg('--cross-mapping-ratio', type=float, metavar='RATIO', default=None,
		help='Ratio for detection of cross-mapping artifacts. Default: %(default)s')
	arg('--min-count', metavar='N', type=int, default=None,
		help='Omit candidates with fewer than N exact occurrences in the input table. '
			'Default: 10 for D, 100 for V and J')

	# --gene=D options
	arg('--d-core-length', metavar='L', type=int, default=6,
		help='Use only D core regions that have at least length L (only '
			'applies when --gene=D). Default: %(default)s')
	arg('--d-core', type=slice_arg, default=slice(2, -2),
		help='D core region location (only applies when --gene=D). '
			'Default: %(default)s')

	arg('--fasta', help='Write discovered sequences to FASTA file')
	arg('table', help='Table with parsed and filtered IgBLAST results')


class Candidate:
	__slots__ = ('name', 'sequence', 'exact_occ', 'max_count', 'other_genes', 'db_name',
		'db_distance', 'cdr3s', 'missing')

	def __init__(self, name, sequence, exact_occ=0, max_count=0, cdr3s=None, other_genes=None,
			db_name=None, db_distance=None):
		self.name = name
		self.sequence = sequence
		self.exact_occ = exact_occ
		self.max_count = max_count  # an upper bound for the exact_occ
		self.cdr3s = cdr3s if cdr3s is not None else set()
		self.other_genes = other_genes if other_genes is not None else set()
		self.db_name = db_name
		self.db_distance = db_distance
		self.missing = ''

	@property
	def unique_CDR3(self):
		return len(self.cdr3s)

	def __repr__(self):
		return 'Candidate({sequence!r}, exact_occ={exact_occ}, max_count={max_count}, ...)'.format(
			sequence=self.sequence,
			exact_occ=self.exact_occ,
			max_count=self.max_count,
		)


class OverlappingSequenceMerger(Merger):
	"""
	Merge sequences that overlap
	"""
	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. If they should not be merged,
		None is returned.
		"""
		m = merge_overlapping(s.sequence, t.sequence)
		if m is not None:
			return Candidate(s.name, m, max_count=t.max_count + s.max_count)

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
				if v.unique_CDR3 == 0:
					continue
				ratio = u.unique_CDR3 / v.unique_CDR3
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
			total_occ = (s.exact_occ + t.exact_occ)
			if total_occ == 0:
				return None
			for u, v in [(s, t), (t, s)]:
				ratio = u.exact_occ / total_occ
				if ratio < self._cross_mapping_ratio:
					# u is probably a cross-mapping artifact of the higher-expressed v
					logger.info('%r is a cross-mapping artifact of %r (ratio %.4f)',
						u.name, v.name, ratio)
					return v

		return None


def count_occurrences(candidates, table_path, search_columns, other_gene, other_errors, merge):
	"""
	Count how often each candidate sequence occurs in the input table.
	The columns named in search_columns are concatenated and searched.
	This circumvents inaccurate IgBLAST alignment boundaries.
	Rows where the column named by other_errors is not zero are ignored.

	The input table is read in chunks to reduce memory usage.

	candidates -- list of candidates
	table_path -- path to input table
	search_columns -- which columns to search. The contained strings are
	    concatenated and then searched.
	merge -- If True, stop searching for other candidates in a single row
	    after one candidate has been found.

	The following attributes of the candidates are updated:

	- exact_occ
	- other_genes
	- cdr3s

	Return the updated list of candidates.
	"""
	candidates_map = {c.sequence: c for c in candidates}

	# Speed up search by looking for the most common sequences first
	search_order = sorted(candidates_map, key=lambda s: candidates_map[s].max_count, reverse=True)
	cols = [other_gene, 'V_errors', 'J_errors', 'CDR3_nt'] + search_columns
	for chunk in pd.read_csv(table_path, usecols=cols, chunksize=10000, sep='\t'):
		fix_columns(chunk)
		chunk = chunk[chunk[other_errors] == 0]
		# concatenate search columns
		if len(chunk) == 0:  # TODO that this is needed is possibly a pandas bug
			continue
		chunk['haystack'] = chunk.loc[:, search_columns].astype(str).sum(axis=1)
		chunk['haystack'] = chunk['haystack'].str.replace('(', '').replace(')', '')

		for row in chunk.itertuples():
			for needle in search_order:
				if needle in row.haystack:
					candidate = candidates_map[needle]
					candidate.exact_occ += 1  # TODO += row.count?
					candidate.other_genes.add(getattr(row, other_gene))
					candidate.cdr3s.add(row.CDR3_nt)
					if merge:
						# When overlapping candidates have been merged,
						# there will be no other pattern that is a
						# substring of the current search pattern.
						break
	return candidates_map.values()


def discard_substring_occurrences(seq_occ_pairs):
	# Sort by sequence length
	seq_occ_pairs = sorted(seq_occ_pairs, key=lambda r: len(r[0]))
	result = []
	while seq_occ_pairs:
		seq, exact_occ = seq_occ_pairs.pop()
		# Keep only those candidates that are not substrings of seq
		tmp = []
		for s, n in seq_occ_pairs:
			if seq.find(s) != -1:
				exact_occ += n
			else:
				tmp.append((s, n))
		seq_occ_pairs = tmp
		result.append((seq, exact_occ))
	return result


def sequence_candidates(table, column, minimum_length, core=slice(None, None)):
	"""
	Generate candidates by clustering all sequences in a column
	(usually V_nt, J_nt). At least two occurrences are required.

	core -- a slice object. If given, the strings in the column are
	       sliced before being clustered.
	"""
	# for sequence, occ in table[column].str[core].value_counts().items():
	# 	if len(sequence) >= minimum_length and occ >= 2:
	# 		yield Candidate(None, sequence, max_count=occ)
	candidates = []
	for sequence, group in table.groupby(column):
		if len(sequence) < minimum_length or len(group) < 2:
			continue
		candidates.append((sequence, len(group)))
	for seq, count in discard_substring_occurrences(candidates):
		yield Candidate(None, seq, max_count=count)


def print_table(candidates, other_gene, missing):
	columns = ['name', 'exact_occ', other_gene + 's', 'CDR3s', 'database', 'database_diff', 'sequence']
	if missing:
		columns.append('missing')
	print(*columns, sep='\t')
	for candidate in candidates:
		columns = [
			candidate.name,
			candidate.exact_occ,
			len(set(candidate.other_genes)),
			candidate.unique_CDR3,
			candidate.db_name if candidate.db_name is not None else '',
			candidate.db_distance if candidate.db_distance is not None else -1,
			candidate.sequence
		]
		if missing:
			columns.append(candidate.missing)
		print(*columns, sep='\t')


def main(args):
	if args.database:
		with FastaReader(args.database) as fr:
			database = list(fr)
		logger.info('Read %d sequences from %r', len(database), args.database)
	else:
		database = None
	column = {'V': 'V_nt', 'J': 'J_nt', 'D': 'D_region'}[args.gene]
	other = 'V' if args.gene in ('D', 'J') else 'J'
	other_gene = other + '_gene'
	other_errors = other + '_errors'
	table = read_table(args.table,
		usecols=['count', 'V_gene', 'J_gene', 'V_errors', 'J_errors', column, 'CDR3_nt'])
	logger.info('Table with %s rows read', len(table))
	table = table[table[other_errors] == 0]
	logger.info('Keeping %s rows that have no %s mismatches', len(table), other)

	if args.merge is None:
		args.merge = args.gene == 'D'
	if args.min_count is None:
		args.min_count = 10 if args.gene == 'D' else 100

	if args.gene == 'D':
		candidates = sequence_candidates(
			table, column, minimum_length=args.d_core_length, core=args.d_core)
	else:
		candidates = sequence_candidates(
			table, column, minimum_length=MINIMUM_CANDIDATE_LENGTH)
	candidates = list(candidates)
	logger.info('Collected %s unique %s sequences', len(candidates), args.gene)

	if args.merge:
		logger.info('Merging overlapping sequences ...')
		# Merge candidate sequences that overlap. If one candidate is longer than
		# another, this is typically a sign that IgBLAST has not extended the
		# alignment long enough.
		merger = OverlappingSequenceMerger()
		for candidate in candidates:
			merger.add(candidate)
		logger.info('After merging overlapping %s sequences, %s remain', args.gene, len(merger))
		candidates = list(merger)
	del table

	logger.info('Counting occurrences ...')

	if args.gene == 'D':
		search_columns = ['VD_junction', 'D_region', 'DJ_junction']
	else:
		search_columns = ['genomic_sequence']
	records = count_occurrences(candidates, args.table, search_columns, other_gene, other_errors,
		args.merge)

	logger.info('%d records', len(records))
	# Assign names etc.
	if database:
		for record in records:
			distances = [(edit_distance(db.sequence, record.sequence), db) for db in database]
			record.db_distance, closest = min(distances, key=lambda x: x[0])
			record.db_name = closest.name

			if record.db_distance == 0:
				record.name = closest.name
			else:
				# Exact db sequence not found, is there one that contains
				# this candidate as a substring?
				is_substring = [db for db in database
					if db.sequence.find(record.sequence) != -1]
				if is_substring:
					db = is_substring[-1]
					start = db.sequence.find(record.sequence)
					prefix = db.sequence[:start]
					suffix = db.sequence[start + len(record.sequence):]
					record.name = db.name
					record.missing = '{}...{}'.format(prefix, suffix)
				else:
					record.name = unique_name(closest.name, record.sequence)
	else:
		for record in records:
			record.name = unique_name(args.gene, record.sequence)

	# Filter by allele ratio
	if args.allele_ratio or args.cross_mapping_ratio:
		arm = AlleleRatioMerger(args.allele_ratio, args.cross_mapping_ratio)
		arm.extend(records)
		records = list(arm)
		logger.info('After filtering by allele ratio and/or cross-mapping ratio, %d candidates remain',
			len(records))

	records = sorted(records, key=lambda r: (r.exact_occ, r.sequence), reverse=True)
	records = [r for r in records if r.exact_occ >= args.min_count]
	print_table(records, other_gene, missing=args.gene == 'D')

	if args.fasta:
		with open(args.fasta, 'w') as f:
			for record in sorted(records, key=lambda r: r.name):
				print('>{}\n{}'.format(record.name, record.sequence), file=f)

	logger.info('Wrote %d genes', len(records))
