"""
Cluster upstream sequences (UTR and leader) for each gene

For each gene, look at all the sequences assigned to it. Take the upstream
sequences and compute a consensus for them. Only those assigned sequences are
taken that have a very low error rate for the V gene match.

Output a FASTA file that contains one consensus sequences for each gene.
"""
import logging
from sqt.align import multialign, consensus
from .table import read_table
#from .utils import iterative_consensus

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('upstream', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=upstream_command)
	subparser.add_argument('--max-error-percentage', '-e', metavar='PERCENT',
		type=float, default=1,
		help='allow PERCENT errors. Default: %(default)s%%.')
	subparser.add_argument('--consensus-threshold', '-t', metavar='PERCENT',
		type=float, default=80,
		help='Threshold for consensus computation. Default: %(default)s%%.')
	subparser.add_argument('--part', choices=['UTR', 'leader', 'UTR+leader'],
		default='UTR+leader', help='Which part of the sequence before the V '
		'gene match to analyze. Default: %(default)s')

	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	return subparser


def upstream_command(args):
	table = read_table(args.table)
	logger.info('%s rows read', len(table))
	table = table[table.V_SHM <= args.max_error_percentage]
	logger.info('%s rows remain after discarding V%%SHM > %s%%', len(table), args.max_error_percentage)
	table['UTR+leader'] = table['UTR'] + table['leader']
	table = table[table[args.part] != '']

	#table = table.loc[:,('name', 'V_gene', 'J_gene', 'V_nt', 'CDR3_nt', 'V_SHM', 'J_SHM', 'UTR', 'leader')].copy()

	for name, group in table.groupby('V_gene'):
		sequences = group[args.part]
		if len(sequences) == 0:
			logger.info('Gene %s has %s assignments, but no usable sequences, skipping.', name, len(group))
			continue
		maxlen = max(len(s) for s in sequences)
		sequences = sequences[[(len(s) >= 0.9 * maxlen) for s in sequences]]
		if len(sequences) == 1:
			cons = sequences.iloc[0]
		else:
			# Keep only those sequences whose length is at least 90% of the longest one
			aligned = multialign(sequences, program='muscle-medium')
			cons = consensus(aligned, threshold=0.8)
		logger.info('Gene %s has %s assignments, %s usable (%s unique sequences). '
			'Consensus has %s N bases.', name, len(group), len(sequences),
			len(set(sequences)), cons.count('N'))
		print('>{} {}_consensus\n{}'.format(name, args.part, cons))
