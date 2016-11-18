"""
Filter V gene candidates (germline and pre-germline filter)

After candidates for novel V genes have been found with the 'discover'
subcommand, this script is used to filter the candidates and make sure that
only true germline genes remain ("germline filter" and "pre-germline filter").
The following filtering and processing steps are performed:

* Discard sequences with N bases
* Discard sequences that come from a consensus over too few source sequences 
* Discard sequences with too few unique CDR3s (CDR3s_exact column)
* Discard sequences with too few unique Js (Js_exact column)
* Discard sequences identical to one of the database sequences (if DB given)
* Discard sequences that do not match a set of known good motifs
* Discard sequences that contain a stop codon (has_stop column)
* Discard near-duplicate sequences

If you provide a whitelist of sequences, then the candidates that appear on it
* are not checked for the cluster size criterion,
* do not need to match a set of known good motifs,
* are never considered near-duplicates,
* are allowed to contain a stop codon.

The filtered table is written to standard output.
"""
import logging
from collections import namedtuple
import pandas as pd
from sqt import FastaReader
from sqt.align import edit_distance

from .utils import UniqueNamer, Merger

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--cluster-size', type=int, metavar='N', default=0,
		help='Consensus must represent at least N sequences. '
		'Default: %(default)s')
	arg('--max-differences', type=int, metavar='MAXDIFF', default=1,
		help='Merge sequences if they have at most MAXDIFF differences. '
		'The one with more CDR3s is kept. Default: %(default)s')
	arg('--cross-mapping-ratio', type=float, metavar='RATIO', default=0.02,
		help='Ratio for detection of cross-mapping artifacts. Default: %(default)s')
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
	arg('--looks-like-V', action='store_true', default=False,
		help='Sequences must look like V genes (uses the looks_like_V column). '
		'Default: Column is ignored')
	arg('--allow-stop', action='store_true', default=False,
		help='Allow stop codons in sequences (uses the has_stop column).'
			'Default: Do not allow stop codons.')
	arg('--whitelist', metavar='FASTA', default=[], action='append',
		help='Sequences that are never discarded or merged with others, '
			'even if criteria for discarding them would apply (except cross-mapping artifact '
			'removal, which is always performed).')
	arg('--fasta', metavar='FILE', help='Write new database in FASTA format to FILE')
	arg('tables', metavar='CANDIDATES.TAB',
		help='Tables (one or more) created by the "discover" command',
		nargs='+')


SequenceInfo = namedtuple('SequenceInfo', 'sequence name CDR3s_exact cluster_size whitelisted is_database cluster_size_is_accurate row')


class SequenceMerger(Merger):
	"""
	Merge sequences that are sufficiently similar into single entries.
	"""
	def __init__(self, max_differences, cross_mapping_ratio):
		super().__init__()
		self._max_differences = max_differences
		self._cross_mapping_ratio = cross_mapping_ratio

	def merged(self, s, t):
		"""
		Given two SequenceInfo objects, decide whether to discard one of them and which one.
		This is used for merging similar candidate sequences.

		Two sequences are considered to be similar if their edit distance is at most
		max_differences (see constructor). If one of the sequences is longer, the 'overhanging'
		bases are ignored at either the 5' end or the 3' end, whichever gives the lower
		edit distance.

		Sequences that are whitelisted are never discarded.
		If two sequences are similar and none of them is whitelisted, then the one with the
		higher number of unique CDR3s is kept.

		Return None if both objects should be kept.
		Return the object to keep otherwise.
		"""
		if len(s.sequence) > len(t.sequence):
			s, t = t, s  # make s always the shorter sequence
		s_seq = s.sequence
		if len(s_seq) != len(t.sequence):
			t_prefix = t.sequence[:len(s_seq)]
			t_suffix = t.sequence[-len(s_seq):]
			dist_prefix = edit_distance(s_seq, t_prefix, max(self._max_differences, 1))
			dist_suffix = edit_distance(s_seq, t_suffix, max(self._max_differences, 1))
			dist = min(dist_prefix, dist_suffix)
		else:
			dist = edit_distance(s_seq, t.sequence, max(self._max_differences, 1))

		# Check for possible cross-mapping
		if self._cross_mapping_ratio and dist == 1 and s.is_database and t.is_database:
			total_count = (s.cluster_size + t.cluster_size)
			for u, v in [(s, t), (t, s)]:
				if u.cluster_size_is_accurate and u.cluster_size / total_count < self._cross_mapping_ratio:
					# u is probably a cross-mapping artifact of the higher-expressed v
					return v

		if dist > self._max_differences:
			return None  # keep both

		if s.whitelisted and t.whitelisted:
			return None  # keep both
		if s.whitelisted:
			return s
		if t.whitelisted:
			return t

		# No sequence is whitelisted if we arrive here
		if s.CDR3s_exact >= t.CDR3s_exact:
			return s
		if len(s_seq) < len(t.sequence):
			return t
		return t


def main(args):
	merger = SequenceMerger(args.max_differences, args.cross_mapping_ratio)

	whitelist = dict()
	for path in args.whitelist:
		for record in FastaReader(path):
			whitelist[record.sequence.upper()] = record.name
	logger.info('%d unique sequences in whitelist', len(whitelist))

	def whitelist_dist(sequence):
		if sequence in whitelist:
			return 0, whitelist[sequence]
		mindist = len(sequence)
		distances = []
		for seq, name in whitelist.items():
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

	# Read in tables
	total_unfiltered = 0
	overall_table = None
	for path in args.tables:
		table = pd.read_csv(path, sep='\t')
		# TODO remove this after deprecation period
		table.rename(columns=dict(
			consensus_seqs='cluster_size',
			window_seqs='cluster_size',
			subset_seqs='cluster_size'), inplace=True)

		i = list(table.columns).index('consensus')
		# whitelist_diff distinguishes between 0 and !=0 only
		# at this point. Accurate edit distances are computed later.
		whitelist_diff = [ (0 if s in whitelist else -1) for s in table['consensus'] ]
		table.insert(i, 'whitelist_diff', pd.Series(whitelist_diff, index=table.index, dtype=int))
		table.insert(i+1, 'closest_whitelist', pd.Series('', index=table.index))

		unfiltered_length = len(table)
		table = table[table.database_diff >= args.minimum_db_diff]
		if 'N_bases' in table.columns:
			table = table[table.N_bases <= args.maximum_N]
		table = table[table.CDR3s_exact >= args.unique_CDR3]
		table = table[table.Js_exact >= args.unique_J]
		if args.looks_like_V:
			table = table[(table.looks_like_V == 1) | (table.whitelist_diff == 0)]
		if not args.allow_stop:
			table = table[(table.has_stop == 0) | (table.whitelist_diff == 0)]
		table = table[(table.cluster_size >= args.cluster_size) | (table.whitelist_diff == 0)]
		table = table.dropna()
		logger.info('Table read from %r contains %s candidate V gene sequences. '
			'%s remain after filtering', path,
			unfiltered_length, len(table))
		if args.whitelist:
			logger.info('Of those, %d are protected by the whitelist', sum(table.whitelist_diff == 0))
		total_unfiltered += unfiltered_length

		if overall_table is None:
			overall_table = table
		else:
			overall_table.append(table)
	del table
	if len(args.tables) > 1:
		logger.info('Read %s tables with %s entries total. '
			'After filtering, %s entries remain.', len(args.tables),
			total_unfiltered, len(overall_table))

	def cluster_size_is_accurate(row):
		return bool(set(row.cluster.split(';')) & set(['all', 'db']))

	for _, row in overall_table.iterrows():
		merger.add(SequenceInfo(row['consensus'], row['name'], row['CDR3s_exact'],
			row['cluster_size'], row['whitelist_diff'] == 0, row['database_diff'] == 0,
			cluster_size_is_accurate(row),
			row.name))  # row.name is the index of the row. It is not row['name'].

	# Discard near-duplicates
	overall_table['is_duplicate'] = pd.Series(True, index=overall_table.index, dtype=bool)
	for info in merger:
		overall_table.loc[info.row, 'is_duplicate'] = False
	overall_table = overall_table[~overall_table.is_duplicate].copy()
	del overall_table['is_duplicate']

	# Name sequences
	overall_table['name'] = overall_table['name'].apply(UniqueNamer())
	overall_table.sort_values(['name'], inplace=True)

	# Because whitelist_dist() is expensive, this is run when
	# all of the filtering has already been done
	if whitelist:
		for row in overall_table.itertuples():
			distance, name = whitelist_dist(overall_table.loc[row[0], 'consensus'])
			overall_table.loc[row[0], 'closest_whitelist'] = name
			overall_table.loc[row[0], 'whitelist_diff'] = distance
	else:
		overall_table.whitelist_diff.replace(-1, '', inplace=True)
	print(overall_table.to_csv(sep='\t', index=False, float_format='%.1f'), end='')

	if args.fasta:
		with open(args.fasta, 'w') as f:
			for _, row in overall_table.iterrows():
				print('>{}\n{}'.format(row['name'], row['consensus']), file=f)

	logger.info('%d sequences in new database', len(overall_table))
