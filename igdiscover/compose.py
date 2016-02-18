"""
Filter V gene candidates (germline/pre-germline filter)

After candidates for novel V genes have been found with the 'discover'
subcommand, this script is used to filter the candidates and make sure that
only true germline genes remain ("germline filter" and "pre-germline filter").
The following filtering and processing steps are performed:

* Discard sequences with N bases
* Discard sequences that come from a consensus over too few source sequences 
* Discard sequences with too few unique CDR3s (CDR3s_exact column)
* Discard sequences identical to one of the database sequences (if DB given)
* Discard sequences that do not match a set of known good motifs
* Merge nearly identical sequences (allowing length differences) into single entries
"""
import logging
import re
from collections import namedtuple
from itertools import zip_longest
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
	arg('--minimum-db-diff', '-b', type=int, metavar='N', default=0,
		help='Sequences must have at least N differences to the database '
		'sequence. Default: %(default)s')
	arg('--maximum-N', '-N', type=int, metavar='COUNT', default=0,
		help='Sequences must have at most COUNT "N" bases. Default: %(default)s')
	arg('--unique-CDR3', '--CDR3s', type=int, metavar='N', default=1,
		help='Sequences must have at least N CDR3s within exact sequence matches. '
		'Default: %(default)s')
	arg('--looks-like-V', action='store_true', default=False,
		help='Sequences must look like V genes (uses the looks_like_V column). '
		'Default: Column is ignored')
	arg('--database', metavar='DATABASE.FASTA',
		help='Existing (to be augmented) database in FASTA format')
	arg('tables', metavar='DISCOVER.TAB',
		help='Tables (one or more) created by the "discover" command',
		nargs='+')


SequenceInfo = namedtuple('SequenceInfo', 'sequence name CDR3s_exact')


class SequenceMerger(Merger):
	"""
	Merge sequences where one is a prefix of the other into single entries.
	"""
	def __init__(self, max_differences):
		super().__init__()
		self._max_differences = max_differences

	def merged(self, s, t):
		"""
		Merge two sequences if one is the prefix of the other. If they should
		not be merged, None is returned.

		s and t must have attributes sequence and name.
		"""
		s_seq = s.sequence
		t_seq = t.sequence
		# Make both sequences the same length - cheap trick to not penalize
		# end gaps
		s_seq += t_seq[len(s_seq):]
		t_seq += s_seq[len(t_seq):]
		if edit_distance(s_seq, t_seq) <= self._max_differences:
			if s.CDR3s_exact >= t.CDR3s_exact:
				return s
			return t
		return None


def main(args):
	merger = SequenceMerger(args.max_differences)
	previous_n = 0
	if args.database:
		for record in FastaReader(args.database):
			previous_n += 1
			merger.add(SequenceInfo(record.sequence.upper(), record.name, 0))  # TODO zero?

	# Read in tables
	total_unfiltered = 0
	tables = []
	for path in args.tables:
		table = pd.read_csv(path, sep='\t')
		# TODO remove this after deprecation period
		table.rename(columns=dict(consensus_seqs='cluster_size', window_seqs='cluster_size', subset_seqs='cluster_size'), inplace=True)

		unfiltered_length = len(table)
		table = table[table.database_diff >= args.minimum_db_diff]
		if 'N_bases' in table.columns:
			table = table[table.N_bases <= args.maximum_N]
		table = table[table.CDR3s_exact >= args.unique_CDR3]
		if args.looks_like_V:
			table = table[table.looks_like_V == 1]
		table = table[table.cluster_size >= args.cluster_size]
		table = table.dropna()
		logger.info('Table read from %r contains %s candidate V gene sequences. '
			'%s remain after filtering', path,
			unfiltered_length, len(table))
		total_unfiltered += unfiltered_length
		tables.append(table)
	if len(args.tables) > 1:
		logger.info('Read %s tables with %s entries total. '
			'After filtering, %s entries remain.', len(args.tables),
			total_unfiltered, sum(map(len, tables)))

	for table in tables:
		for _, row in table.iterrows():
			merger.add(SequenceInfo(row['consensus'], row['name'], row['CDR3s_exact']))

	namer = UniqueNamer()
	n = 0
	for info in merger:
		n += 1
		print('>{}\n{}'.format(namer(info.name), info.sequence))

	if args.database:
		logger.info('Old database had %s sequences, new database has %s sequences (difference: %s)', previous_n, n, n - previous_n)
	else:
		logger.info('New database has %s sequences', n)
