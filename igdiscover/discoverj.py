"""
Rudimentary discovery of J genes

Prints out the most frequent J sequences in FASTA format
"""
import logging
from collections import Counter

from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('table', help='Table with parsed and filtered IgBLAST results')


def main(args):
	table = read_table(args.table)
	table = table[table.V_errors == 0]
	sequences = Counter(table.J_nt)
	for i, (seq, freq) in enumerate(sequences.most_common(), 1):
		if freq < 10:
			break
		print('>J_{} freq={}\n{}'.format(i, freq, seq))
