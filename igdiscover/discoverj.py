"""
Rudimentary discovery of J genes

Print out the most frequent J sequences in FASTA format.
If one sequence is suffix of another, the longer version is kept.
"""
import logging
from collections import Counter, namedtuple
import pandas as pd
from .utils import Merger
from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('table', help='Table with parsed and filtered IgBLAST results')


SequenceInfo = namedtuple('SequenceInfo', 'name sequence groups frequencies frequency')


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
        # make s the longer sequence
        if len(s.sequence) < len(t.sequence):
            s, t = t, s
        if s.sequence.endswith(t.sequence): # or s.sequence.startswith(t.sequence):

            frequencies = s.frequencies.copy()
            lendiff = len(s.sequence) - len(t.sequence)
            assert lendiff > 0, (s, t)
            for index, freq in t.frequencies.items():
                assert index + lendiff not in frequencies
                frequencies[index + lendiff] = freq
            return SequenceInfo('name', s.sequence, s.groups + t.groups, frequencies, t.frequency + s.frequency)
        else:
            return None


def main(args):
	table = read_table(args.table)
	table = table[table.V_errors == 0]

	sequences = Counter(table.J_nt)
	logger.info('Table contains %s rows with %s distinct J sequences', len(table), len(sequences))

	merger = SequenceMerger()
	# Note that the merging result depends on the order in which we iterate.
	# groupby returns keys in sorted order.
	for sequence, group in table.groupby('J_nt'):
		merger.add(SequenceInfo('name', sequence, [group], {0: len(group)}, len(group)))
	logger.info('After merging, %s sequences remain', len(merger))

	print('sequence', 'count', 'V_genes', 'CDR3s', 'name', sep='\t')
	i = 0
	for record in sorted(merger, key=lambda r: r.sequence):
		group = pd.concat(record.groups)
		if record.frequency < 100:  # TODO
			continue
		names = Counter(group.J_gene)
		if len(names) > 1:
			name = ', '.join('{}({})'.format(name, names[name]) for name in sorted(names))
		else:
			name = next(iter(names))
		print(record.sequence,
			sum(record.frequencies.values()),
			len(set(group.V_gene)),
			len(set(group.CDR3_nt)),
			name,
			sep='\t')
		i += 1
		if i == 100:
			break
