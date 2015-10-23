"""
Cluster upstream sequences (UTR and leader/before the V gene match)
"""
import logging

from .table import read_table
#from .utils import iterative_consensus

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('discover', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--threads', '-j', type=int, default=min(4, available_cpu_count()),
		help='Number of threads. Default: no. of available CPUs, but at most 4')
	subparser.add_argument('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors. Default: %(default)s.')
	subparser.add_argument('--consensus-threshold', '-t', metavar='PERCENT', type=float, default=60,
		help='Threshold for consensus computation. Default: %(default)s%%.')

	...

	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser




def upstream_command(args):
	table = read_table(args.table)
	table = table.loc[:,('name', 'V_gene', 'J_gene', 'V_nt', 'CDR3_nt', 'V_SHM', 'J_SHM')].copy()


	TO DO

	d['UTR+leader'] = d['UTR'] + d['leader']

	max_error_percentage = 1
	for part in ('leader', 'UTR', 'UTR+leader'):
		print()
		print('##### Analyzing', part)
		print()
		for name, group in d.groupby('V_gene'):
			sequences = group[(group.V_SHM <= max_error_percentage) & (group[part] != '')][part]
			if len(sequences) == 0:
				print('Gene', name, 'has', len(group), 'assignments, but no usable sequences, skipping.')
				continue
			maxlen = max(len(s) for s in sequences)
			sequences = sequences[[(len(s) >= 0.9 * maxlen) for s in sequences]]
			#print('min, max:', min(len(s) for s in sequences), max(len(s) for s in sequences))
			if len(sequences) == 1:
				cons = sequences[0]
			else:
				# Keep only those sequences whose length is at least 90% of the longest one
				aligned = multialign(sequences, program='muscle-medium')
				cons = consensus(aligned, threshold=0.8)
			print('Gene', name, 'has', len(group), 'assignments.', len(sequences), 'usable sequences.', len(set(sequences)), 'unique sequences. Consensus:')
			print(cons)


	# Discard rows with any mutation within J at all
	logger.info('%s rows read', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	logger.info('Using an error rate window of %.1f%% to %.1f%%', args.left, args.right)
	logger.info('Approximate comparisons between V gene sequence and consensus allow %.1f%% errors.', v_error_rate*100)




	logger.info('%s consensus sequences for %s gene(s) computed', n_consensus, len(groups))








