"""
Cluster upstream sequences (UTR and leader) for each gene

For each gene, look at all the sequences assigned to it. Take the upstream
sequences and compute a consensus for them. Only those assigned sequences are
taken that have a very low error rate for the V gene match.

Output a FASTA file that contains one consensus sequence for each gene.
"""
import logging
from sqt.align import multialign, consensus
from .table import read_table
from collections import Counter
#from .utils import iterative_consensus

logger = logging.getLogger(__name__)

# When computing a UTR consensus, ignore sequences that deviate by more than
# this factor from the median length.
UTR_MEDIAN_DEVIATION = 0.1


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('upstream', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=upstream_command)
	subparser.add_argument('--max-error-percentage', '-e', metavar='PERCENT',
		type=float, default=1,
		help='allow PERCENT errors. Default: %(default)s%%.')
	subparser.add_argument('--consensus-threshold', '-t', metavar='PERCENT',
		type=float, default=75,
		help='Threshold for consensus computation. Default: %(default)s%%.')
	subparser.add_argument('--part', choices=['UTR', 'leader', 'UTR+leader'],
		default='UTR+leader', help='Which part of the sequence before the V '
		'gene match to analyze. Default: %(default)s')
	subparser.add_argument('--debug', default=False, action='store_true',
		help='Enable debugging output')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	return subparser


def upstream_command(args):
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
		if len(group) == 0:
			logger.info('Gene %s has %s assignments, but no usable sequences, skipping.', name, len(group))
			continue
		median_length = group['UTR_length'].median()

		#counter = Counter(group['UTR_length'])
		## If there is anything at freq. 10 or higher, take the longest of those
		#frequent = [length for length, count in counter.items() if count >= 10]
		#if frequent:
			#longest = max(frequent)
		#else:

		lower = median_length * (1.0 - UTR_MEDIAN_DEVIATION)
		upper = median_length * (1.0 + UTR_MEDIAN_DEVIATION)
		counter = Counter(len(s) for s in group['UTR'])
		logger.debug('Sequence lengths (length: count): %s',  ', '.join('{}: {}'.format(l,c) for l, c in counter.most_common()))
		logger.debug('Median: %s. Lower bound: %s. Upper bound: %s', median_length, lower, upper)
		sequences = group[args.part]
		sequences = sequences[group['UTR_length'] >= lower]
		sequences = sequences[group['UTR_length'] <= upper]
		if len(sequences) == 0:
			logger.info('Gene %s has %s assignments, but lengths are too different, skipping.', name, len(group))
			continue
		assert len(sequences) > 0
		if len(sequences) == 1:
			cons = sequences.iloc[0]
		else:
			# Keep only those sequences whose length is at least 90% of the longest one
			aligned = multialign(sequences, program='muscle-medium')
			cons = consensus(aligned, threshold=args.consensus_threshold/100)
		logger.info('Gene %s has %s assignments, %s usable (%s unique sequences). '
			'Consensus has %s N bases.', name, len(group), len(sequences),
			len(set(sequences)), cons.count('N'))
		if cons.count('N') > 0:
			n_consensus_with_n += 1
		n_written += 1
		print('>{} {}_consensus\n{}'.format(name, args.part, cons))

	logger.info('Wrote a consensus for %s of %s genes (%s had N bases)', n_written, n, n_consensus_with_n)
