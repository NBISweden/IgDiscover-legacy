"""
Discover new V genes within a single antibody library.
"""
import logging
import os.path
from collections import Counter, OrderedDict
import numpy as np
from sqt.align import multialign, consensus, edit_distance
from .table import read_table

logger = logging.getLogger(__name__)

MINGROUPSIZE_CONSENSUS = 10


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors (default: %(default)s)')
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Use "all" to compute for all genes.')
	subparser.add_argument('--left', '-l', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at least this %%SHM (default: %(default)s)', default=0)
	subparser.add_argument('--right', '-r', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at most this %%SHM (default: %(default)s)', default=100)
	subparser.add_argument('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser


def sister_sequence(group, program='muscle-medium'):
	"""

	"""
	sequences = OrderedDict()
	# TODO Perhaps create the dict in such a way that those with the most
	# abundant no. of errors are put in first.
	for _, row in group.iterrows():
		sequences[row.name] = row.V_nt
	aligned = multialign(sequences, program=program)
	cons = consensus(aligned, threshold=0.6)
	return cons.strip('N')


def discover_command(args):
	v_error_rate = args.error_rate / 100
	assert 0 <= v_error_rate <= 1
	table = read_table(args.table)
	#table = table.loc[:,['V_gene', 'V%SHM', 'V_nt', 'name']]
	#t = table.loc[:,('name', 'V_gene', 'V_nt', 'V%SHM')]

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	genes = set(args.gene)
	n = 0
	logger.info('Using a %%SHM window of %.1f%% to %.1f%%', args.left, args.right)
	logger.info('Approximate comparisons between V gene sequence and consensus allow %.1f%% errors.', v_error_rate*100)
	for gene, group in table.groupby('V_gene'):
		if not ('all' in genes or gene in genes):
			continue
		if len(group) < MINGROUPSIZE_CONSENSUS:
			logger.info('Skipping %s as the number of sequences is too small (%s)', gene, len(group))
			continue
		logger.info('Working on gene %s', gene)

		# Create various subsets of the full group
		group_in_shm_range = group[(group.V_SHM >= args.left) & (group.V_SHM <= args.right)]
		s = sister_sequence(group_in_shm_range)

		group = group.copy()
		group['consensus_diff'] = [ edit_distance(v_nt, s) for v_nt in group.V_nt ]
		group_exact_V = group[group.V_nt == s]
		group_approximate_V = group[group.consensus_diff <= len(s) * v_error_rate]

		for description, g in (
				('sequences in total were assigned to this gene', group),
				('sequences are within the %SHM window (consensus was computed from these)', group_in_shm_range),
				('sequences match the consensus exactly', group_exact_V),
				('sequences match the consensus approximately', group_approximate_V)):
			logger.info('%s %s (%.1f%%):', len(g), description, len(g) / len(group) * 100)
			logger.info('   %s unique J genes used', len(set(g.J_gene)))
			logger.info('   %s unique CDR3 sequences used', len(set(g.CDR3_nt)))
		logger.info('Consensus computed from sequences within %%SHM window has %s “N” bases', s.count('N'))

		# Print this last so it doesn’t mess up output too bad in case stdout
		# isn’t redirected anywhere.
		print('>{}_sister\n{}'.format(gene, s))

		if args.table_output and len(group_approximate_V) > 0:
			if not os.path.exists(args.table_output):
				os.mkdir(args.table_output)
			path = os.path.join(args.table_output, gene + '.tab')
			group_approximate_V.sort('consensus_diff').to_csv(path, sep='\t')
			logger.info('Wrote %s', path)

		n += 1
	if genes and n > 1:
		logger.info('%s consensus sequences computed', n)
