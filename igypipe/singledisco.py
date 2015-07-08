"""
Discover new V genes within a single antibody library. Existing V sequences are
grouped by their V gene assignment and within each group, consensus sequences
are computed.
"""
import csv
import logging
import sys
import os.path
from collections import Counter, OrderedDict, namedtuple
import numpy as np
from sqt import SequenceReader
from sqt.align import multialign, consensus, edit_distance
from .table import read_table

logger = logging.getLogger(__name__)

MINGROUPSIZE_CONSENSUS = 10


Groupinfo = namedtuple('Groupinfo', 'count unique_J unique_CDR3')


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors. Default: %(default)s.')
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Default: Compute for all genes.')
	subparser.add_argument('--left', '-l', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at least this error rate (in percent). Default: %(default)s', default=0)
	subparser.add_argument('--right', '-r', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at most this error rate (in percent). Default: %(default)s', default=100)
	subparser.add_argument('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	subparser.add_argument('--database', metavar='FASTA', default=None,
		help='FASTA file with V genes. If provided, differences between consensus and database will be computed.')

	subparser.add_argument('--consensus-output', '-c', metavar='FASTA', default=None,
		help='Output consensus sequences in FASTA format to this file.')
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

	if args.database:
		with SequenceReader(args.database) as sr:
			database = { record.name: record.sequence.upper() for record in sr }
	else:
		database = dict()

	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	logger.info('Using an error rate window of %.1f%% to %.1f%%', args.left, args.right)
	logger.info('Approximate comparisons between V gene sequence and consensus allow %.1f%% errors.', v_error_rate*100)

	if args.consensus_output:
		consensus_output = open(args.consensus_output, 'w')
	else:
		consensus_output = None
	writer = csv.writer(sys.stdout, delimiter='\t')
	writer.writerow([
		'gene',
		'total_seqs',
		'total_unique_J',
		'total_unique_CDR3',
		'window_seqs',
		'window_unique_J',
		'window_unique_CDR3',
		'exact_seqs',
		'exact_unique_J',
		'exact_unique_CDR3',
		'approx_seqs',
		'approx_unique_J',
		'approx_unique_CDR3',
		'N_bases',
		'database_diff',
		'consensus'
	])
	genes = set(args.gene)
	n = 0
	for gene, group in table.groupby('V_gene'):
		if not ('all' in genes or len(genes) == 0 or gene in genes):
			continue
		if len(group) < MINGROUPSIZE_CONSENSUS:
			writer.writerow([gene, len(group)])
			logger.info('Skipping %s as the number of sequences is too small (%s)', gene, len(group))
			continue

		logger.info('Working on gene %s', gene)

		# Create various subsets of the full group
		group_in_shm_range = group[(group.V_SHM >= args.left) & (group.V_SHM <= args.right)]
		sister = sister_sequence(group_in_shm_range)

		group = group.copy()
		group['consensus_diff'] = [ edit_distance(v_nt, sister) for v_nt in group.V_nt ]
		group_exact_V = group[group.V_nt == sister]
		group_approximate_V = group[group.consensus_diff <= len(sister) * v_error_rate]

		info = dict()
		for key, g in (
				('total', group),
				('window', group_in_shm_range),
				('exact', group_exact_V),
				('approx', group_approximate_V)):
			unique_J = len(set(g.J_gene))
			unique_CDR3 = len(set(s for s in g.CDR3_nt if s))
			info[key] = Groupinfo(count=len(g), unique_J=unique_J, unique_CDR3=unique_CDR3)

		if gene in database:
			database_diff = edit_distance(sister, database[gene])
		else:
			database_diff = None
		if consensus_output:
			print('>{}_sister\n{}'.format(gene, sister), file=consensus_output)

		if args.table_output and len(group_approximate_V) > 0:
			if not os.path.exists(args.table_output):
				os.mkdir(args.table_output)
			path = os.path.join(args.table_output, gene + '.tab')
			group_approximate_V.sort('consensus_diff').to_csv(path, sep='\t')
			logger.info('Wrote %s', path)

		row = [gene]
		for key in ('total', 'window', 'exact', 'approx'):
			row.extend([info[key].count, info[key].unique_J, info[key].unique_CDR3])
		row.extend([sister.count('N'), database_diff, sister])
		writer.writerow(row)
		sys.stdout.flush()
		n += 1
	if consensus_output:
		consensus_output.close()
	logger.info('%s consensus sequences computed', n)
