"""
Discover candidate new V genes within a single antibody library.

Existing V sequences are grouped by their V gene assignment, and within each
group, consensus sequences are computed.
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
from .utils import iterative_consensus, sequence_hash, downsampled, SerialPool
from .cluster import cluster_sequences
from .utils import looks_like_V_gene, Merger

logger = logging.getLogger(__name__)

random.seed(123)

MINGROUPSIZE_CONSENSUS = 5
MAXIMUM_SUBSAMPLE_SIZE = 1600

Groupinfo = namedtuple('Groupinfo', 'count unique_J unique_CDR3')

SisterInfo = namedtuple('SisterInfo', 'sequence requested name group')


def add_arguments(parser):
	arg = parser.add_argument
	arg('--threads', '-j', type=int, default=min(4, available_cpu_count()),
		help='Number of threads. Default: no. of available CPUs, but at most 4')
	arg('--consensus-threshold', '-t', metavar='PERCENT', type=float, default=60,
		help='Threshold for consensus computation. Default: %(default)s%%.')
	arg('--prefix', default='', metavar='PREFIX',
		help='Add PREFIX before sequence names')
	arg('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Default: Compute for all genes.')
	arg('--left', '-l', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at least this error rate (in percent). Default: %(default)s', default=0)
	arg('--right', '-r', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at most this error rate (in percent). Default: %(default)s', default=100)
	arg('--window-width', '-w', type=float, metavar='PERCENT',
		help='Compute consensus for all PERCENT-wide windows. Default: do not compute', default=None)
	arg('--cluster', action='store_true', default=False,
		help='Cluster sequences by similarity and compute consensus')
	arg('--max-n-bases', type=int, default=0, metavar='MAXN',
		help='Remove rows that have more than MAXN "N" nucleotides. If >0, an '
			'N_bases column is added. Default: 0')
	arg('--subsample', metavar='N', type=int, default=1000,
		help='When clustering, use N randomly chosen sequences. Default: %(default)s')
	arg('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	arg('--database', metavar='FASTA', default=None,
		help='FASTA file with V genes. If provided, differences between consensus and database will be computed.')
	arg('--consensus-output', '-c', metavar='FASTA', default=None,
		help='Output consensus sequences in FASTA format to this file.')
	arg('--ignore-J', action='store_true', default=False,
		help='Include also rows without J assignment or J%%SHM>0.')
	arg('--approx', action='store_true', default=False,
		help='Count also approximate matches (adds three columns to output table).')
	arg('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors. Default: %(default)s.')
	arg('table', help='Table with parsed IgBLAST results')  # nargs='+'


class SisterMerger(Merger):
	"""
	Merge very similar consensus sequences into single entries. This could be
	seen as a type of clustering using very specific criteria. Two sequences
	are merged if one is the prefix of the other, allowing differences where
	one of the sequences has an 'N' base.
	"""
	def merged(self, s, t):
		chars = []
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
			chars.append(c)
		seq = ''.join(chars)
		requested = s.requested or t.requested
		name = s.name + ';' + t.name
		# take union of groups
		group = pd.concat([s.group, t.group]).groupby(level=0).last()
		return SisterInfo(seq, requested, name, group)


class Discoverer:
	"""
	Discover candidates for novel V genes.
	"""
	def __init__(self, database, windows, left, right, cluster, table_output,
			  prefix, consensus_threshold, v_error_rate, downsample,
			  cluster_subsample_size, approx_columns, max_n_bases):
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
		self.cluster_subsample_size = cluster_subsample_size
		self.approx_columns = approx_columns
		self.max_n_bases = max_n_bases

	def _sister_sequence(self, group):
		"""
		For a given group, compute a consensus sequence over the V gene sequences
		in that group.
		"""
		return iterative_consensus(list(group.V_nt), program='muscle-medium',
			threshold=self.consensus_threshold/100,
			maximum_subsample_size=self.downsample)

	def __call__(self, args):
		"""
		Discover new V genes. args is a tuple (gene, group)
		gene -- name of the gene
		group -- a pandas DataFrame with the group corresponding to the gene
		"""
		gene, group = args
		# Collect all 'sister' sequences (consensus sequences)
		sisters = SisterMerger()
		group = group.copy()
		for left, right in self.windows:
			left, right = float(left), float(right)
			group_in_window = group[(left <= group.V_SHM) & (group.V_SHM < right)]
			if len(group_in_window) < MINGROUPSIZE_CONSENSUS:
				continue
			sister = self._sister_sequence(group_in_window)
			if left == int(left):
				left = int(left)
			if right == int(right):
				right = int(right)
			requested = (left, right) == (self.left, self.right)
			name = '{}-{}'.format(left, right)
			sisters.add(SisterInfo(sister, requested, name, group_in_window))

		if self.cluster:
			indices = downsampled(list(group.index), self.cluster_subsample_size)
			sequences = list(group.V_nt.loc[indices])
			df, linkage, clusters = cluster_sequences(sequences)
			logger.info('Clustering sequences for %r gave %d cluster(s)', gene,
				len(set(clusters)))
			cluster_indices = [ [] for _ in range(max(clusters) + 1) ]
			for i, cluster_id in enumerate(clusters):
				cluster_indices[cluster_id].append(indices[i])

			cl = 1
			for ind in cluster_indices:
				group_in_window = group.loc[ind]
				if len(group_in_window) < MINGROUPSIZE_CONSENSUS:
					continue
				sister = self._sister_sequence(group_in_window)
				name = 'cl{}'.format(cl)
				info = SisterInfo(sister, False, name, group_in_window)
				sisters.add(info)
				cl += 1

		rows = []
		for sister_info in sisters:
			sister = sister_info.sequence
			n_bases = sister.count('N')
			if n_bases > self.max_n_bases:
				continue
			group_exact_V = group[group.V_nt == sister]
			if self.approx_columns:
				group['consensus_diff'] = [ edit_distance(v_nt, sister) for v_nt in group.V_nt ]
				group_approximate_V = group[group.consensus_diff <= len(sister) * self.v_error_rate]

			groups = (
				('window', sister_info.group),
				('exact', group_exact_V))
			if self.approx_columns:
				groups +=  (('approx', group_approximate_V), )
			info = dict()
			for key, g in groups:
				unique_J = len(set(s for s in g.J_gene if s))
				unique_CDR3 = len(set(s for s in g.CDR3_nt if s))
				count = len(g.index)
				info[key] = Groupinfo(count=count, unique_J=unique_J, unique_CDR3=unique_CDR3)
			if gene in self.database:
				database_diff = edit_distance(sister, self.database[gene])
			else:
				database_diff = None

			if self.approx_columns:
				assert info['exact'].count <= info['approx'].count

			# Build the row for the output table
			# TODO use UniqueNamer here
			sequence_id = '{}{}_{}'.format(self.prefix, gene.rsplit('_S', 1)[0], sequence_hash(sister))
			row = [gene, sequence_id, sister_info.name if len(sister_info.group) < len(group) else 'all']
			for key, _ in groups:
				row.extend([info[key].count, info[key].unique_J, info[key].unique_CDR3])
			if self.max_n_bases:
				row += [n_bases]
			row += [database_diff, int(looks_like_V_gene(sister)), sister]
			rows.append(row)

			# If a window was requested via --left/--right, write the 'approx'
			# subset to a separate file.
			if self.table_output and self.approx_columns and any(si.requested for si in sister_info) and len(group_approximate_V) > 0:
				if not os.path.exists(self.table_output):
					os.mkdir(self.table_output)
				path = os.path.join(self.table_output, gene + '.tab')
				group_approximate_V.sort('consensus_diff').to_csv(path, sep='\t')
				logger.info('Wrote %s for window %s-%s', path, self.left, self.right)
		return rows


def main(args):
	v_error_rate = args.error_rate / 100
	assert 0 <= v_error_rate <= 1

	if args.database:
		with SequenceReader(args.database) as sr:
			database = { record.name: record.sequence.upper() for record in sr }
	else:
		database = dict()

	table = read_table(args.table)
	table = table.loc[:,('name', 'V_gene', 'J_gene', 'V_nt', 'CDR3_nt', 'V_SHM', 'J_SHM')].copy()

	logger.info('%s rows read', len(table))
	if not args.ignore_J:
		# Discard rows with any mutation within J at all
		table = table[table.J_SHM == 0][:]
		logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	logger.info('Using an error rate window of %.1f%% to %.1f%%', args.left, args.right)
	if args.approx:
		logger.info('Approximate comparisons between V gene sequence and consensus allow %.1f%% errors.', v_error_rate*100)

	if args.consensus_output:
		consensus_output = open(args.consensus_output, 'w')
	else:
		consensus_output = None
	writer = csv.writer(sys.stdout, delimiter='\t')
	columns = [
		'source_gene',
		'name',
		'consensus_from',
		'consensus_seqs',
		'consensus_unique_J',
		'consensus_unique_CDR3',
		'exact_seqs',
		'exact_unique_J',
		'exact_unique_CDR3',
	]
	if args.approx:
		columns += ['approx_seqs', 'approx_unique_J', 'approx_unique_CDR3']
	if args.max_n_bases:
		columns += ['N_bases']
	columns += [
		'database_diff',
		'looks_like_V',
		'consensus'
	]
	writer.writerow(columns)
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
		args.table_output, args.prefix, args.consensus_threshold, v_error_rate,
		MAXIMUM_SUBSAMPLE_SIZE, cluster_subsample_size=args.subsample,
		approx_columns=args.approx, max_n_bases=args.max_n_bases)
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
