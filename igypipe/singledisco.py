"""
Discover potential new V genes within a single antibody library. Existing V
sequences are grouped by their V gene assignment and within each group,
consensus sequences are computed.
"""
import csv
import logging
import sys
import os.path
from collections import Counter, OrderedDict, namedtuple
import multiprocessing
import numpy as np
import pandas as pd
from sqt import SequenceReader
from sqt.align import multialign, consensus, edit_distance
from sqt.utils import available_cpu_count
from .table import read_table

logger = logging.getLogger(__name__)

MINGROUPSIZE_CONSENSUS = 10

Groupinfo = namedtuple('Groupinfo', 'count unique_J unique_CDR3')


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--threads', '-j', type=int, default=available_cpu_count(),
		help='Number of threads. Default: no. of available CPUs (%(default)s)')
	subparser.add_argument('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors. Default: %(default)s.')
	subparser.add_argument('--consensus-threshold', '-t', metavar='PERCENT', type=float, default=60,
		help='Threshold for consensus computation. Default: %(default)s%%.')
	subparser.add_argument('--prefix', default='', metavar='PREFIX',
		help='Add PREFIX before sequence names')
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Default: Compute for all genes.')
	subparser.add_argument('--left', '-l', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at least this error rate (in percent). Default: %(default)s', default=0)
	subparser.add_argument('--right', '-r', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at most this error rate (in percent). Default: %(default)s', default=100)
	subparser.add_argument('--window-width', '-w', type=float, metavar='PERCENT',
		help='Compute consensus for all PERCENT-wide windows. Default: do not compute', default=None)
	subparser.add_argument('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	subparser.add_argument('--database', metavar='FASTA', default=None,
		help='FASTA file with V genes. If provided, differences between consensus and database will be computed.')
	subparser.add_argument('--consensus-output', '-c', metavar='FASTA', default=None,
		help='Output consensus sequences in FASTA format to this file.')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser


def sister_sequence(group, program='muscle-medium', threshold=0.6):
	"""

	"""
	sequences = OrderedDict()
	# TODO Perhaps create the dict in such a way that those with the most
	# abundant no. of errors are put in first.
	for _, row in group.iterrows():
		sequences[row.name] = row.V_nt
	aligned = multialign(sequences, program=program)
	cons = consensus(aligned, threshold=threshold)
	return cons.strip('N')


def sequence_hash(s):
	"""Return a four-digit hash of a string (S000 to S999)"""
	return 'S{:03}'.format(abs(hash(s)) % 1000)


class Discoverer:
	"""
	Discover candidates for novel V genes.
	"""
	def __init__(self, database, windows, left, right, table_output, prefix, consensus_threshold, v_error_rate):
		self.database = database
		self.windows = windows
		self.left = left
		self.right = right
		self.table_output = table_output
		self.prefix = prefix
		self.consensus_threshold = consensus_threshold
		self.v_error_rate = v_error_rate

	def __call__(self, args):
		gene, group = args
		#logger.info('Working on gene %s', gene)

		# Collect all 'sister' sequences (consensus sequences)
		sisters = OrderedDict()  # sequence -> list of (left, right) tuples
		group = group.copy()
		for left, right in self.windows:
			if left == int(left):
				left = int(left)
			if right == int(right):
				right = int(right)
			group_in_window = group[(left <= group.V_SHM) & (group.V_SHM < right)]
			if len(group_in_window) < MINGROUPSIZE_CONSENSUS:
				continue
			sister = sister_sequence(group_in_window, threshold=self.consensus_threshold/100)
			if sister in sisters:
				sisters[sister].append((left, right, group_in_window))
			else:
				sisters[sister] = [(left, right, group_in_window)]

		rows = []
		for sister, sister_info in sisters.items():
			dists = [ edit_distance(v_nt, sister) for v_nt in group.V_nt ]
			assert len(dists) == len(group)

			group['consensus_diff'] = dists
			group_exact_V = group[group.V_nt == sister]
			group_approximate_V = group[group.consensus_diff <= len(sister) * self.v_error_rate]


			# Instead of concatenating, the groups should be unioned, but not
			# sure how to do this in pandas.
			sister_windows = [(w[0], w[1]) for w in sister_info]
			group_in_window = pd.concat(w[2] for w in sister_info)

			info = dict()
			for key, g in (
					('total', group),  # TODO re-done for every sister
					('window', group_in_window),
					('exact', group_exact_V),
					('approx', group_approximate_V)):
				unique_J = len(set(g.J_gene))
				unique_CDR3 = len(set(s for s in g.CDR3_nt if s))
				# Correct for concatenated dataframes by taking the set of indices.
				count = len(set(g.index))
				info[key] = Groupinfo(count=count, unique_J=unique_J, unique_CDR3=unique_CDR3)
			if gene in self.database:
				database_diff = edit_distance(sister, self.database[gene])
			else:
				database_diff = None
			n_bases = sister.count('N')
			window_str = ';'.join('{}-{}'.format(l, r) for l, r in sister_windows)
			sequence_id = '{}{}_{}'.format(self.prefix, gene, sequence_hash(sister))

			# Build the row for the output table
			row = [gene, window_str]
			for key in ('total', 'window', 'exact', 'approx'):
				row.extend([info[key].count, info[key].unique_J, info[key].unique_CDR3])
			row.extend([n_bases, database_diff, sequence_id, sister])
			rows.append(row)

			# If a window was requested via --left/--right, write the 'approx'
			# subset to a separate file.
			if self.table_output and (self.left, self.right) in sister_windows and len(group_approximate_V) > 0:
				if not os.path.exists(self.table_output):
					os.mkdir(self.table_output)
				path = os.path.join(self.table_output, gene + '.tab')
				group_approximate_V.sort('consensus_diff').to_csv(path, sep='\t')
				logger.info('Wrote %s for window %s-%s', path, self.left, self.right)
		return rows


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
		'window',
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
		'name',
		'consensus'
	])
	genes = set(args.gene)
	if args.window_width:
		windows = [ (start, start + args.window_width) for start in np.arange(0, 20, args.window_width) ]
	else:
		windows = []
	windows.append((args.left, args.right))

	groups = []
	for gene, group in table.groupby('V_gene'):
		if not ('all' in genes or len(genes) == 0 or gene in genes):
			continue
		if len(group) < MINGROUPSIZE_CONSENSUS:
			continue
		groups.append((gene, group))

	discoverer = Discoverer(database, windows, args.left, args.right,
		args.table_output, args.prefix, args.consensus_threshold, v_error_rate)
	n_consensus = 0
	with multiprocessing.Pool(args.threads) as pool:
		for rows in pool.imap(discoverer, groups, chunksize=1):
			writer.writerows(rows)
			sys.stdout.flush()

			if consensus_output:
				for row in rows:
					print('>{} window:{}\n{}'.format(row[-2], row[1], row[-1]), file=consensus_output)
				n_consensus += 1
	if consensus_output:
		consensus_output.close()
	logger.info('%s consensus sequences for %s gene(s) computed', n_consensus, len(groups))
