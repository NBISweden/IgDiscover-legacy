"""
Create new V gene database from V gene candidates.

After potentially new V gene sequences have been discovered with the discover
subcommand, this script can be used to create a new V gene database. The
following filtering and processing steps are performed:

* Discard sequences with N bases
* Discard sequences with too few unique CDR3s (exact_unique_CDR3 column)
* Discard sequences identical to one of the database sequences
* Discard sequences that do not match a set of known good motifs
* Merge nearly identical sequences (allowing length differences) into single entries
"""
import logging
import re
from collections import namedtuple
from itertools import zip_longest
import pandas as pd
from sqt import FastaReader
from sqt.dna import amino_acid_regex
from .utils import UniqueNamer

logger = logging.getLogger(__name__)


def add_arguments(parser):
	parser.add_argument('--minimum-db-diff', '-b', type=int, metavar='DIST', default=0,
		help='Sequences must have at least DIST differences to the database sequence. Default: %(default)s')
	parser.add_argument('--maximum-N', '-N', type=int, metavar='COUNT', default=0,
		help='Sequences must have at most COUNT N bases. Default: %(default)s')
	parser.add_argument('--unique-CDR3', type=int, metavar='COUNT', default=5,
		help='Sequences must have at least COUNT exact unique CDR3s. Default: %(default)s')
	parser.add_argument('--database', metavar='DATABASE.FASTA',
		help='Existing (to be augmented) database in FASTA format')
	parser.add_argument('tables', metavar='DISCOVER.TAB',
		help='Table (zero or more) created by the "discover" command', nargs='*')


SequenceInfo = namedtuple('SequenceInfo', 'sequence name')


class Merger:
	"""
	Merge sequences where one is a prefix of the other into single entries.

	TODO unify with discover.SisterMerger?
	"""
	def __init__(self):
		self._sequences = []

	def add(self, info):
		# See if we already have a similar sequence
		for i, s in enumerate(self._sequences):
			c = self._merged(s, info)
			if c is not None:
				self._sequences[i] = c
				break
		else:
			self._sequences.append(info)

	def __iter__(self):
		yield from self._sequences

	@staticmethod
	def _merged(s, t):
		"""
		Merge two sequences if one is the prefix of the other. If they should
		not be merged, None is returned.

		s and t must have attributes sequence and name.
		"""
		seq = []
		for c1, c2 in zip_longest(s.sequence, t.sequence):
			if c1 is None:
				c = c2
			elif c2 is None:
				c = c1
			#elif c1 == 'N':
				#c = c2
			#elif c2 == 'N':
				#c = c1
			elif c1 != c2:
				return None
			else:
				assert c1 == c2
				c = c1
			seq.append(c)
		seq = ''.join(seq)
		name = ';'.join(set(s.name.split(';')).union(t.name.split(';')))
		return SequenceInfo(seq, name)


def main(args):
	merger = Merger()
	previous_n = 0
	if args.database:
		for record in FastaReader(args.database):
			previous_n += 1
			merger.add(SequenceInfo(record.sequence.upper(), record.name))

	# Read in tables
	total_unfiltered = 0
	tables = []
	for path in args.tables:
		table = pd.read_csv(path, sep='\t')
		unfiltered_length = len(table)
		table = table[table.database_diff >= args.minimum_db_diff]
		if 'N_bases' in table.columns:
			table = table[table.N_bases <= args.maximum_N]
		table = table[table.exact_unique_CDR3 >= args.unique_CDR3]
		table = table[[looks_like_V_gene(s) for s in table.consensus]]

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
			merger.add(SequenceInfo(row['consensus'], row['name']))

	namer = UniqueNamer()
	n = 0
	for info in merger:
		n += 1
		print('>{}\n{}'.format(namer(info.name), info.sequence))

	if args.database:
		logger.info('Old database had %s sequences, new database has %s sequences (difference: %s)', previous_n, n, n - previous_n)
	else:
		logger.info('New database has %s sequences', n)
