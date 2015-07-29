"""
Discover potential new V genes within a single antibody library. Existing V
sequences are grouped by their V gene assignment and within each group,
consensus sequences are computed.
"""
import csv
import logging
import sys
import os.path
from collections import OrderedDict, namedtuple
import multiprocessing
import random
from itertools import zip_longest

import numpy as np
import pandas as pd
from sqt import SequenceReader
from sqt.align import edit_distance
from sqt.utils import available_cpu_count
from .table import read_table
from .utils import iterative_consensus, sequence_hash, downsampled
from .cluster import cluster_sequences, cluster_consensus

logger = logging.getLogger(__name__)

random.seed(123)

MINGROUPSIZE_CONSENSUS = 10
MAXIMUM_SUBSAMPLE_SIZE = 1600
CLUSTER_SUBSAMPLE_SIZE = 300

Groupinfo = namedtuple('Groupinfo', 'count unique_J unique_CDR3')

SisterInfo = namedtuple('SisterInfo', 'sequence requested name group')


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--threads', '-j', type=int, default=min(4, available_cpu_count()),
		help='Number of threads. Default: no. of available CPUs, but at most 4')
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
	subparser.add_argument('--cluster', action='store_true', default=False,
		help='Cluster sequences by similarity and compute consensus')
	subparser.add_argument('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	subparser.add_argument('--database', metavar='FASTA', default=None,
		help='FASTA file with V genes. If provided, differences between consensus and database will be computed.')
	subparser.add_argument('--consensus-output', '-c', metavar='FASTA', default=None,
		help='Output consensus sequences in FASTA format to this file.')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser


class SerialPool:
	"""
	An alternative to multiprocessing.Pool that runs things in parallel for
	easier debugging
	"""
	def __init__(self, *args, **kwargs):
		pass

	def __enter__(self):
		return self

	def __exit__(self, *args):
		pass

	def imap(self, func, iterable, chunksize):
		for i in iterable:
			yield func(i)


def sister_sequence(group, program='muscle-medium', threshold=0.6, maximum_subsample_size=1600):
	"""
	For a given group, compute a consensus sequence over the V gene sequences
	in that group.
	"""
	return iterative_consensus(list(group.V_nt), program, threshold, maximum_subsample_size=maximum_subsample_size)


class SisterMerger:
	"""
	Merge very similar consensus sequences into single entries. This could be
	seen as a type of clustering using very specific criteria. Two sequences
	are merged if one is the prefix of the other, allowing differences where
	one of the sequences has an 'N' base.
	"""
	def __init__(self):
		self.sisters = []

	def add(self, info):
		self.sisters.append(info)

	def _merge_all(self):
		if not self.sisters:
			return []
		merged = [self.sisters[0]]
		for s in self.sisters[1:]:
			for i, m in enumerate(merged):
				c = self._merged(m, s)
				if c is not None:
					merged[i] = c
					break
			else:
				# Found no similar sister sequence
				merged.append(s)
		return merged

	def __iter__(self):
		for m in self._merge_all():
			yield m

	@staticmethod
	def _merged(s, t):
		seq = []
		for c1, c2 in zip_longest(s.sequence, t.sequence):
			if c1 is None:
				c = c2
			elif c2 is None:
				c = c1
			elif c1 == 'N':
				c = c2
			elif c2 == 'N':
				c = c1
			elif c1 != c2:
				return None
			else:
				assert c1 == c2
				c = c1
			seq.append(c)
		seq = ''.join(seq)
		requested = s.requested or t.requested
		name = s.name + ';' + t.name
		# take union of groups
		group = pd.concat([s.group, t.group]).groupby(level=0).last()
		return SisterInfo(seq, requested, name, group)


class Discoverer:
	"""
	Discover candidates for novel V genes.
	"""
	def __init__(self, database, windows, left, right, cluster, table_output, prefix, consensus_threshold, v_error_rate, downsample):
		self.database = database
		self.windows = windows
		self.left = left
		self.right = right
		self.cluster = cluster
		self.table_output = table_output
		self.prefix = prefix
		self.consensus_threshold = consensus_threshold
		self.v_error_rate = v_error_rate
		self.downsample = downsample

	def __call__(self, args):
		gene, group = args
		# Collect all 'sister' sequences (consensus sequences)
		sisters = SisterMerger()
		group = group.copy()
		for left, right in self.windows:
			left, right = float(left), float(right)
			group_in_window = group[(left <= group.V_SHM) & (group.V_SHM < right)]
			if len(group_in_window) < MINGROUPSIZE_CONSENSUS:
				continue
			sister = sister_sequence(group_in_window, threshold=self.consensus_threshold/100, maximum_subsample_size=self.downsample)
			if left == int(left):
				left = int(left)
			if right == int(right):
				right = int(right)
			requested = (left, right) == (self.left, self.right)
			name = '{}-{}'.format(left, right)
			sisters.add(SisterInfo(sister, requested, name, group_in_window))

		if self.cluster:
			indices = downsampled(list(group.index), CLUSTER_SUBSAMPLE_SIZE)
			sequences = list(group.V_nt.loc[indices])
			df, linkage, clusters = cluster_sequences(sequences)
			for i, sister in enumerate(cluster_consensus(sequences, clusters), 1):
				name = 'cl{}'.format(i)
				info = SisterInfo(sister, False, name, group.loc[indices])
				sisters.add(info)

		rows = []
		for sister_info in sisters:
			sister = sister_info.sequence
			dists = [ edit_distance(v_nt, sister) for v_nt in group.V_nt ]
			assert len(dists) == len(group)

			group['consensus_diff'] = dists
			group_exact_V = group[group.V_nt == sister]
			group_approximate_V = group[group.consensus_diff <= len(sister) * self.v_error_rate]

			info = dict()
			for key, g in (
					('total', group),  # TODO re-done for every sister
					('window', sister_info.group),
					('exact', group_exact_V),
					('approx', group_approximate_V)):
				unique_J = len(set(g.J_gene))
				unique_CDR3 = len(set(s for s in g.CDR3_nt if s))
				# Correct for concatenated dataframes by taking the set of indices.
				# TODO
				assert len(set(g.index)) == len(g.index)
				count = len(set(g.index))
				info[key] = Groupinfo(count=count, unique_J=unique_J, unique_CDR3=unique_CDR3)
			if gene in self.database:
				database_diff = edit_distance(sister, self.database[gene])
			else:
				database_diff = None
			n_bases = sister.count('N')

			# Build the row for the output table
			sequence_id = '{}{}_{}'.format(self.prefix, gene, sequence_hash(sister))
			row = [gene, sister_info.name]
			for key in ('total', 'window', 'exact', 'approx'):
				row.extend([info[key].count, info[key].unique_J, info[key].unique_CDR3])
			row.extend([n_bases, database_diff, sequence_id, sister])
			rows.append(row)

			# If a window was requested via --left/--right, write the 'approx'
			# subset to a separate file.
			if self.table_output and any(si.requested for si in sister_info) and len(group_approximate_V) > 0:
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
	table = table.loc[:,('name', 'V_gene', 'J_gene', 'V_nt', 'CDR3_nt', 'V_SHM', 'J_SHM')].copy()

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

	discoverer = Discoverer(database, windows, args.left, args.right, args.cluster,
		args.table_output, args.prefix, args.consensus_threshold, v_error_rate, MAXIMUM_SUBSAMPLE_SIZE)
	n_consensus = 0

	Pool = SerialPool if args.threads == 1 else multiprocessing.Pool
	with Pool(args.threads) as pool:
		for rows in pool.imap(discoverer, groups, chunksize=1):
			writer.writerows(rows)
			sys.stdout.flush()

			if consensus_output:
				for row in rows:
					print('>{} window:{}\n{}'.format(row[-2], row[1], row[-1]), file=consensus_output)
			n_consensus += len(rows)
	if consensus_output:
		consensus_output.close()
	logger.info('%s consensus sequences for %s gene(s) computed', n_consensus, len(groups))
