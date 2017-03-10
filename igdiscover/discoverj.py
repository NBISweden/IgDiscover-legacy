"""
Rudimentary discovery of J genes

Print out the most frequent J sequences in FASTA format.
If one sequence is suffix of another, the longer version is kept.
"""
import logging
from collections import Counter, namedtuple
from .utils import Merger
from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('table', help='Table with parsed and filtered IgBLAST results')


SequenceInfo = namedtuple('SequenceInfo', 'name sequence frequency')


class SequenceMerger(Merger):
	"""
	Merge sequences where one is a suffix of the other.
	"""
	def __init__(self):
		super().__init__()

	def merged(self, s, t):
		"""
		Merge two sequences if one is a suffix of the other. If they should
		not be merged, None is returned.

		s and t must have attributes sequence and frequency
		"""
		if s.sequence.endswith(t.sequence):
			return SequenceInfo('name', s.sequence, s.frequency + t.frequency)
		elif t.sequence.endswith(s.sequence):
			return SequenceInfo('name', t.sequence, s.frequency + t.frequency)
		else:
			return None


def main(args):
	table = read_table(args.table)
	table = table[table.V_errors == 0]
	merger = SequenceMerger()
	sequences = Counter(table.J_nt)
	for sequence, frequency in sequences.items():
		merger.add(SequenceInfo('name', sequence, frequency))

	i = 1
	for record in merger:
		if record.frequency < 10:
			continue
		print('>J_{} freq={}\n{}'.format(i, record.frequency, record.sequence))
		i += 1
