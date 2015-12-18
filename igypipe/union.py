"""
Compute union of sequences in multiple FASTA files
"""
import logging
from collections import namedtuple
from sqt import FastaReader
from sqt.align import edit_distance

from .utils import Merger

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--max-differences', type=int, metavar='MAXDIFF', default=0,
		help='Merge sequences if they have at most MAXDIFF differences. '
		' Default: %(default)s')
	arg('fasta', help='FASTA file', nargs='+')


SequenceInfo = namedtuple('SequenceInfo', 'sequence name')


class SequenceMerger(Merger):
	"""
	Merge sequences where one is a prefix of the other into single entries.
	"""
	def __init__(self, max_differences):
		super().__init__()

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
		if s_seq == t_seq:
			return s
		return None

def main(args):
	merger = SequenceMerger(args.max_differences)
	n_read = 0
	for path in args.fasta:
		n = 0
		for record in FastaReader(path):
			merger.add(SequenceInfo(record.sequence.upper(), record.name))
			n += 1
		n_read += n
		logger.info('Read %s sequences from %s', n, path)
	logger.info('Read %s sequences from %s files', n_read, len(args.fasta))

	for info in merger:
		print('>{}\n{}'.format(info.name, info.sequence))
