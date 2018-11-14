"""
Filter V gene candidates (germline and pre-germline filter)

After candidates for novel V genes have been found with the 'discover'
subcommand, this script is used to filter the candidates and make sure that
only true germline genes remain ("germline filter" and "pre-germline filter").
The following filtering and processing steps are performed on each candidate
separately:

* Discard sequences with N bases
* Discard sequences that come from a consensus over too few source sequences (unless whitelisted)
* Discard sequences with too few unique CDR3s (CDR3_clusters column)
* Discard sequences with too few unique Js (Js_exact column)
* Discard sequences identical to one of the database sequences (if DB given)
* Discard sequences that contain a stop codon (has_stop column) (unless whitelisted)

The following criteria involve comparison of candidates against each other:

* Discard sequences that are too similar to another (unless whitelisted)
* Discard sequences that are cross-mapping artifacts
* Discard sequences that have a too low allele ratio
* Discard sequences with too few unique Ds relative to other alleles (Ds_exact column)

If you provide a whitelist of sequences, then the candidates that appear on it
* are not checked for the cluster size criterion,
* are never considered near-duplicates,
* are allowed to contain a stop codon.

The filtered table is written to standard output.
"""
import sys
import logging
from collections import namedtuple
import pandas as pd
from sqt import FastaReader
from sqt.align import edit_distance

from .utils import UniqueNamer, Merger, is_same_gene, ChimeraFinder

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--cluster-size', type=int, metavar='N', default=0,
		help='Consensus must represent at least N sequences. '
		'Default: %(default)s')
	arg('--cross-mapping-ratio', type=float, metavar='RATIO', default=0.02,
		help='Ratio for detection of cross-mapping artifacts. Default: %(default)s')
	arg('--clonotype-ratio', '--allele-ratio', type=float, metavar='RATIO', default=0.1,
		help='Required ratio of "clonotypes" counts between alleles. '
			'Works only for genes named "NAME*ALLELE". Default: %(default)s')
	arg('--exact-ratio', type=float, metavar='RATIO', default=0.1,
		help='Required ratio of "exact" counts between alleles. '
			'Works only for genes named "NAME*ALLELE". Default: %(default)s')
	arg('--cdr3-shared-ratio', type=float, metavar='RATIO', default=1.0,
		help='Maximum allowed CDR3_shared_ratio. Default: %(default)s')
	arg('--minimum-db-diff', '-b', type=int, metavar='N', default=0,
		help='Sequences must have at least N differences to the database '
		'sequence. Default: %(default)s')
	arg('--maximum-N', '-N', type=int, metavar='COUNT', default=0,
		help='Sequences must have at most COUNT "N" bases. Default: %(default)s')
	arg('--unique-CDR3', '--CDR3s', type=int, metavar='N', default=1,
		help='Sequences must have at least N unique CDR3s within exact sequence matches. '
		'Default: %(default)s')
	# The default for unique-J is 0 because we might work on data without
	# any assigned J genes.
	arg('--unique-J', type=int, metavar='N', default=0,
		help='Sequences must have at least N unique Js within exact sequence matches. '
		'Default: %(default)s')
	arg('--unique-D-ratio', type=float, metavar='RATIO', default=None,
		help='Discard a sequence if another allele of this gene exists '
		'such that the ratio between their Ds_exact is less than RATIO')
	arg('--unique-D-threshold', type=int, metavar='THRESHOLD', default=10,
		help='Apply the --unique-D-ratio filter only if the Ds_exact of the other '
		'allele is at least THRESHOLD')
	arg('--allow-stop', action='store_true', default=False,
		help='Allow stop codons in sequences (uses the has_stop column).'
			'Default: Do not allow stop codons.')
	# arg('--allow-chimeras', action='store_true', default=False,
	# 	help='Do not filter out chimeric sequences. Default: Discard chimeras')
	arg('--whitelist', metavar='FASTA', default=[], action='append',
		help='Sequences that are never discarded or merged with others, '
			'even if criteria for discarding them would apply (except cross-mapping artifact '
			'removal, which is always performed).')
	arg('--fasta', metavar='FILE', help='Write new database in FASTA format to FILE')
	arg('--annotate', metavar='FILE',
		help='Write candidates.tab with filter annotations to FILE')
	arg('tables', metavar='CANDIDATES.TAB',
		help='Tables (one or more) created by the "discover" command',
		nargs='+')


SequenceInfo = namedtuple('_SequenceInfo', ['sequence', 'name', 'clonotypes', 'exact', 'Ds_exact',
	'cluster_size', 'whitelisted', 'is_database', 'cluster_size_is_accurate', 'CDR3_start', 'row'])


class CandidateFilter:
	"""
	Merge sequences that are sufficiently similar into single entries.
	"""
	def __init__(self, cross_mapping_ratio, clonotype_ratio, exact_ratio,
			unique_d_ratio,
			unique_d_threshold):
		self._items = []
		self._cross_mapping_ratio = cross_mapping_ratio
		self._clonotype_ratio = clonotype_ratio
		self._exact_ratio = exact_ratio
		self._unique_d_ratio = unique_d_ratio
		self._unique_d_threshold = unique_d_threshold

	def add(self, item):
		# This method could possibly be made simpler if the graph structure
		# was made explicit.

		items = []
		for existing_item in self._items:
			m = self.merged(existing_item, item)
			if m is None:
				items.append(existing_item)
			else:
				item = m
		items.append(item)
		self._items = items

	def extend(self, iterable):
		for i in iterable:
			self.add(i)

	def __iter__(self):
		if self._items and hasattr(self._items, 'name'):
			yield from sorted(self._items, key=lambda x: x.name)
		else:
			yield from self._items

	def __len__(self):
		return len(self._items)

	def merged(self, s: SequenceInfo, t: SequenceInfo):
		"""
		Given two SequenceInfo objects, decide whether to discard one of them and which one.
		This is used for merging similar candidate sequences.

		Two sequences are Sequences can also be discarded if llele ratio, cross-mapping ratio or unique_d ratio criteria

		If one of the sequences is longer, the 'overhanging'
		bases are ignored at either the 5' end or the 3' end, whichever gives the lower
		edit distance.


		If two sequences are similar and
		Sequences that are whitelisted are not discarded if the only reason for discarding
		them is their similarity
		If two sequences are similar and none of them is whitelisted, then the one with the
		higher number of unique CDR3s is kept.

		Return None if both objects should be kept.
		Return the object to keep otherwise.
		"""
		if len(s.sequence) > len(t.sequence):
			s, t = t, s  # make s always the shorter sequence

		# When computing edit distance between the two sequences, ignore the
		# bases in the 3' end that correspond to the CDR3
		s_no_cdr3 = s.sequence[:s.CDR3_start]
		t_no_cdr3 = t.sequence[:t.CDR3_start]
		if len(s_no_cdr3) != len(t_no_cdr3):
			t_prefix = t_no_cdr3[:len(s_no_cdr3)]
			t_suffix = t_no_cdr3[-len(s_no_cdr3):]
			dist_prefix = edit_distance(s_no_cdr3, t_prefix, 1)
			dist_suffix = edit_distance(s_no_cdr3, t_suffix, 1)
			dist = min(dist_prefix, dist_suffix)
		else:
			dist = edit_distance(s_no_cdr3, t_no_cdr3, 1)

		# Check for cross-mapping
		# We check for dist <= 1 (and not ==1) since there may be
		# differences in the CDR3 parts of the sequences that we
		# ignore when computing edit distance.
		if self._cross_mapping_ratio and dist <= 1 and s.is_database and t.is_database:
			total_count = (s.cluster_size + t.cluster_size)
			for u, v in [(s, t), (t, s)]:
				ratio = u.cluster_size / total_count
				if u.cluster_size_is_accurate and ratio < self._cross_mapping_ratio:
					# u is probably a cross-mapping artifact of the higher-expressed v
					logger.info('%r is a cross-mapping artifact of %r (ratio %.4f)',
						u.name, v.name, ratio)
					return v

		# Check criteria based on whether the two sequences are alleles of the same gene
		if is_same_gene(s.name, t.name):
			# Check clonotype allele ratio. Somewhat similar to cross-mapping, but
			# this uses sequence names to decide whether two genes can be
			# alleles of each other and the ratio is between the
			# CDR3_clusters values
			if self._clonotype_ratio:
				for u, v in [(s, t), (t, s)]:
					if v.clonotypes == 0:
						continue
					ratio = u.clonotypes / v.clonotypes
					if ratio < self._clonotype_ratio:
						logger.info('Clonotype allele ratio %.4f too low for %r compared to %r',
							ratio, u.name, v.name)
						return v

			# Check exact V sequence occurrence allele ratio
			if self._exact_ratio:
				for u, v in [(s, t), (t, s)]:
					if v.exact == 0:
						continue
					ratio = u.exact / v.exact
					if ratio < self._exact_ratio:
						logger.info('Allele ratio of exact occurrences %.4f too low for %r compared to %r',
							ratio, u.name, v.name)
						return v

			if self._unique_d_ratio:
				for u, v in [(s, t), (t, s)]:
					if u.cluster_size >= v.cluster_size and u.Ds_exact >= self._unique_d_threshold:
						ratio = v.Ds_exact / u.Ds_exact
						if ratio < self._unique_d_ratio:
							logger.info('Ds_exact ratio %.4f too low for %r compared to %r',
								ratio, v.name, u.name)
							return u

		if dist > 0:
			return None  # keep both

		# If we get here, then the two candidates are too similar. Unless they are
		# whitelisted, one needs to be discarded.
		if s.whitelisted and t.whitelisted:
			return None  # keep both
		if s.whitelisted:
			return s
		if t.whitelisted:
			return t

		# No sequence is whitelisted if we arrive here
		if s.clonotypes >= t.clonotypes:
			return s
		if len(s.sequence) < len(t.sequence):
			return t
		return t


class Whitelist:
	def __init__(self):
		self._sequences = dict()

	def add_fasta(self, path):
		with FastaReader(path) as fr:
			for record in fr:
				self._sequences[record.sequence.upper()] = record.name

	def closest(self, sequence):
		"""
		Search for the whitelist sequence that is closest to the given sequence.

		Return tuple (distance, name).
		"""
		if sequence in self._sequences:
			return 0, self._sequences[sequence]
		mindist = len(sequence)
		distances = []
		for seq, name in self._sequences.items():
			ed = edit_distance(seq, sequence, maxdiff=mindist)
			distances.append((ed, name))
			if ed == 1:
				# We know ed does not get smaller because the
				# 'sequence in whitelist' check
				# above covers that
				return ed, name
			mindist = min(mindist, ed)
		distance, name = min(distances)
		return distance, name

	def __len__(self):
		return len(self._sequences)

	def __contains__(self, other):
		return other in self._sequences


def is_chimera(table, whitelist):
	result = pd.Series('', index=table.index, dtype=object)
	if whitelist:
		whitelisted = table[table['whitelist_diff'] == 0]
	else:
		whitelisted = table[table['database_diff'] == 0]
	chimera_finder = ChimeraFinder(list(whitelisted['consensus']))
	for row in table[(table['whitelist_diff'] != 0) & (table['database_diff'] != 0)].itertuples():
		query = row.consensus
		chimera_result = chimera_finder.find_exact(query)
		if chimera_result:
			(prefix_length, prefix_indices, suffix_indices) = chimera_result
			suffix_length = len(query) - prefix_length
			prefix_name = whitelisted.iloc[prefix_indices[0]]['name']
			suffix_row = whitelisted.iloc[suffix_indices[0]]
			suffix_name = suffix_row['name']
			suffix_sequence = suffix_row['consensus']
			# logger.info('Candidate %s (diffs.: %d) appears to be a chimera of '
			# 	'%s:1..%d and %s:%d..%d', row.name, row.whitelist_diff, prefix_name, prefix_length,
			# 	suffix_name, len(suffix_sequence) - suffix_length + 1, len(suffix_sequence))
			result.loc[row.Index] = '{}:1..{}+{}:{}..{}'.format(prefix_name, prefix_length,
				suffix_name, len(suffix_sequence) - suffix_length + 1, len(suffix_sequence))

	return result


def mark_rows(table, condition, reason):
	table.loc[condition, 'why_filtered'] += reason + ';'
	table.loc[condition, 'is_filtered'] += 1
	logger.info('Marked %s candidates as %r', sum(condition), reason)


def main(args):
	if args.unique_D_threshold <= 1:
		sys.exit('--unique-D-threshold must be at least 1')
	merger = CandidateFilter(
		cross_mapping_ratio=args.cross_mapping_ratio,
		clonotype_ratio=args.clonotype_ratio,
		exact_ratio=args.exact_ratio,
		unique_d_ratio=args.unique_D_ratio,
		unique_d_threshold=args.unique_D_threshold
	)
	whitelist = Whitelist()
	for path in args.whitelist:
		whitelist.add_fasta(path)
	logger.info('%d unique sequences in whitelist', len(whitelist))

	# Read in tables
	total = 0
	overall_table = None
	for path in args.tables:
		table = pd.read_csv(path, sep='\t')
		i = list(table.columns).index('consensus')
		# whitelist_diff distinguishes between 0 and !=0 only
		# at this point. Accurate edit distances are computed later.
		whitelist_diff = [(0 if s in whitelist else -1) for s in table['consensus']]
		# TODO rename to is_whitelisted
		table.insert(i, 'whitelist_diff', pd.Series(whitelist_diff, index=table.index, dtype=int))
		table.insert(i+1, 'closest_whitelist', pd.Series('', index=table.index))
		table.insert(3, 'why_filtered', pd.Series('', index=table.index))
		table.insert(3, 'is_filtered', pd.Series(0, index=table.index))

		mark_rows(table, table.database_diff < args.minimum_db_diff, 'dbdiff_low')

		if 'N_bases' in table.columns:
			mark_rows(table, table.N_bases > args.maximum_N, 'too_many_N_bases')
		mark_rows(table, table.CDR3s_exact < args.unique_CDR3, 'too_low_CDR3s_exact')
		mark_rows(table, table.CDR3_shared_ratio > args.cdr3_shared_ratio, 'too_high_CDR3_shared_ratio')
		mark_rows(table, table.Js_exact < args.unique_J, 'too_low_Js_exact')
		if not args.allow_stop:
			mark_rows(table, (table.has_stop != 0) & (table.whitelist_diff != 0), 'has_stop')
		mark_rows(table, (table.cluster_size < args.cluster_size) & (table.whitelist_diff != 0),
			'too_low_cluster_size')
		table['database_changes'].fillna('', inplace=True)
		# table = table.dropna()  TODO why was this in here?
		logger.info('Table read from %r contains %s candidate V gene sequences. '
			'%s remain after per-entry filtering', path,
			len(table), sum(table.is_filtered == 0))
		if args.whitelist:
			logger.info('Of those, %d are protected by the whitelist', sum(table.whitelist_diff == 0))
		total += len(table)

		if overall_table is None:
			overall_table = table
		else:
			overall_table.append(table)
	del table
	if len(args.tables) > 1:
		logger.info('Read %s tables with %s entries total. '
			'After per-entry filtering, %s entries remain.', len(args.tables),
			len(overall_table), sum(overall_table.is_filtered == 0))

	def cluster_size_is_accurate(row):
		return bool(set(row.cluster.split(';')) & {'all', 'db'})

	for _, row in overall_table.iterrows():
		if row['is_filtered'] > 0:
			continue
		merger.add(SequenceInfo(
			sequence=row['consensus'],
			name=row['name'],
			clonotypes=row['clonotypes'],
			exact=row['exact'],
			Ds_exact=row['Ds_exact'],
			cluster_size=row['cluster_size'],
			whitelisted=row['whitelist_diff'] == 0,
			is_database=row['database_diff'] == 0,
			cluster_size_is_accurate=cluster_size_is_accurate(row),
			CDR3_start=row.get('CDR3_start', 10000),  # TODO backwards compatibility
			row=row.name,  # row.name is the index of the row. It is not row['name'].
		))

	# Discard near-duplicates
	overall_table.loc[overall_table.is_filtered > 0, 'is_duplicate'] = False
	overall_table.loc[overall_table.is_filtered == 0, 'is_duplicate'] = True
	for info in merger:
		overall_table.loc[info.row, 'is_duplicate'] = False
	mark_rows(overall_table, overall_table.is_duplicate, 'is_duplicate')
	del overall_table['is_duplicate']

	# Name sequences
	overall_table['name'] = overall_table['name'].apply(UniqueNamer())
	overall_table.sort_values(['name'], inplace=True)

	# Because whitelist_dist() is expensive, this is run when
	# all of the filtering has already been done
	if whitelist:
		for row in overall_table.itertuples():
			# TODO skipping this is just a performance optimization
			if row.is_filtered > 0:
				continue
			distance, name = whitelist.closest(overall_table.loc[row[0], 'consensus'])
			overall_table.loc[row[0], 'closest_whitelist'] = name
			overall_table.loc[row[0], 'whitelist_diff'] = distance
	else:
		overall_table.whitelist_diff.replace(-1, '', inplace=True)

	i = list(overall_table.columns).index('database_diff')
	# TODO
	# chimeras are computed for the full table, not the filtered one. what should be done?
	overall_table.insert(i, 'chimera', is_chimera(overall_table, whitelist))

	# Discard chimeric sequences
	# if not args.allow_chimeras:
	#   overall_table = overall_table[~is_chimera(overall_table, whitelist)].copy()

	filtered_table = overall_table[overall_table.is_filtered == 0]
	del filtered_table['is_filtered']
	del filtered_table['why_filtered']
	print(filtered_table.to_csv(sep='\t', index=False, float_format='%.2f'), end='')

	if args.annotate:
		with open(args.annotate, 'w') as f:
			print(overall_table.to_csv(sep='\t', index=False, float_format='%.2f'), end='', file=f)
	if args.fasta:
		with open(args.fasta, 'w') as f:
			for _, row in filtered_table.iterrows():
				print('>{}\n{}'.format(row['name'], row['consensus']), file=f)

	logger.info('%d sequences in new database', len(filtered_table))
