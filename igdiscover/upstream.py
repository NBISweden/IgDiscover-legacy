"""
Cluster upstream sequences (UTR and leader) for each gene

For each gene, look at all the sequences assigned to it. Take the upstream
sequences and compute a consensus for them. Only those assigned sequences are
taken that have a very low error rate for the V gene match.

Output a FASTA file that contains one consensus sequence for each gene.
"""
import logging
from collections import Counter
from sqt.align import multialign, consensus

from .table import read_table
from .utils import iterative_consensus

logger = logging.getLogger(__name__)

# When computing a UTR consensus, ignore sequences that deviate by more than
# this factor from the median length.
UTR_MEDIAN_DEVIATION = 0.1

def add_arguments(parser):
	arg = parser.add_argument
	arg('--max-error-percentage', '-e', metavar='PERCENT',
		type=float, default=1,
		help='Allow PERCENT errors in V gene match. Default: %(default)s%%.')
	arg('--consensus-threshold', '-t', metavar='PERCENT',
		type=float, default=75,
		help='Threshold for consensus computation. Default: %(default)s%%.')
	arg('--part', choices=['UTR', 'leader', 'UTR+leader'],
		default='UTR+leader', help='Which part of the sequence before the V '
		'gene match to analyze. Default: %(default)s')
	arg('--debug', default=False, action='store_true',
		help='Enable debugging output')
	arg('table', help='Table with parsed IgBLAST results')


def main(args):
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	table = read_table(args.table)
	logger.info('%s rows read', len(table))
	table = table[table.V_SHM <= args.max_error_percentage]
	logger.info('%s rows remain after discarding V%%SHM > %s%%', len(table), args.max_error_percentage)
	table['UTR+leader'] = table['UTR'] + table['leader']
	table = table[table[args.part] != '']
	table['UTR_length'] = [ len(s) for s in table['UTR'] ]

	n = 0
	n_written = 0
	n_consensus_with_n = 0
	for name, group in table.groupby('V_gene'):
		n += 1
		assert len(group) != 0
		counter = Counter(group['UTR_length'])
		logger.debug('Sequence lengths (length: count): %s',  ', '.join('{}: {}'.format(l,c) for l, c in counter.most_common()))

		if args.part == 'leader':
			sequences = list(group[args.part])
		else:
			# Take all sequences that are at least as long as the tenth longest
			# sequence. This gives us at least ten usable sequences, enough for
			# a consensus.
			length_threshold = sorted(group['UTR_length'], reverse=True)[:10][-1]
			sequences = list(group[group['UTR_length'] >= length_threshold][args.part])

		if len(sequences) == 0:
			logger.info('Gene %s has %s assignments, but lengths are too different, skipping.', name, len(group))
			continue
		assert len(sequences) > 0
		if len(sequences) == 1:
			cons = sequences[0]
		else:
			# Keep only those sequences whose length is at least 90% of the longest one
			cons = iterative_consensus(sequences, program='muscle-medium', threshold=args.consensus_threshold/100)
			# If the sequence lengths are different, the beginning can be unclear
			cons = cons.lstrip('N')
		logger.info('Gene %s has %s assignments, %s usable (%s unique sequences). '
			'Consensus has %s N bases.', name, len(group), len(sequences),
			len(set(sequences)), cons.count('N'))
		if cons.count('N') > 0:
			n_consensus_with_n += 1
		n_written += 1
		print('>{} {}_consensus\n{}'.format(name, args.part, cons))

	logger.info('Wrote a consensus for %s of %s genes (%s had N bases)', n_written, n, n_consensus_with_n)
