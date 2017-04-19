"""
Discover candidate new V genes within a single antibody library.

Existing V sequences are grouped by their V gene assignment, and within each
group, consensus sequences are computed.
"""
import csv
import hashlib
import logging
import sys
import os.path
from collections import OrderedDict, namedtuple, Counter
import multiprocessing
import random
from itertools import zip_longest
import numpy as np
import pandas as pd
from sqt import SequenceReader
from sqt.align import edit_distance
from sqt.utils import available_cpu_count

from .table import read_table
from .utils import iterative_consensus, unique_name, downsampled, SerialPool, Merger, has_stop
from .cluster import cluster_sequences
from .species import looks_like_V_gene

logger = logging.getLogger(__name__)


MINGROUPSIZE = 5
MAXIMUM_SUBSAMPLE_SIZE = 1600


def add_arguments(parser):
	arg = parser.add_argument
	arg('--threads', '-j', type=int, default=min(4, available_cpu_count()),
		help='Number of threads. Default: no. of available CPUs, but at most 4')
	arg('--seed', type=int, default=None,
		help='Seed value for random numbers for reproducible runs.')
	arg('--consensus-threshold', '-t', metavar='PERCENT', type=float, default=60,
		help='Threshold for consensus computation. Default: %(default)s%%')
	arg('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. '
			'Default: Compute for all genes.')
	arg('--left', '-l', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at least this '
			'error rate (in percent). Default: %(default)s', default=0)
	arg('--right', '-r', type=float, metavar='ERROR-RATE',
		help='For consensus, include only sequences that have at most this '
			'error rate (in percent). Default: %(default)s', default=100)
	arg('--window-width', '-w', type=float, metavar='PERCENT',
		help='Compute consensus for all PERCENT-wide windows. Set to 0 to '
			'disable. Default: %(default)s', default=2)
	arg('--no-cluster', dest='cluster', action='store_false', default=True,
		help='Do not run linkage cluster analysis.')
	arg('--max-n-bases', type=int, default=0, metavar='MAXN',
		help='Remove rows that have more than MAXN "N" nucleotides. If >0, an '
			'N_bases column is added. Default: %(default)s')
	arg('--subsample', metavar='N', type=int, default=1000,
		help='When clustering, use N randomly chosen sequences. Default: %(default)s')
	arg('--table-output', '-o', metavar='DIRECTORY',
		help='Output tables for all analyzed genes to DIRECTORY. '
			'Files will be named <GENE>.tab.')
	arg('--database', metavar='FASTA', default=None,
		help='FASTA file with V genes. If provided, differences between consensus '
			'and database will be computed.')
	arg('--ignore-J', action='store_true', default=False,
		help='Include also rows without J assignment or J%%SHM>0.')
	arg('--approx', action='store_true', default=False,
		help='Count also approximate matches (adds three columns to output table).')
	arg('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors. '
			'Default: %(default)s')
	arg('--exact-copies', metavar='N', type=int, default=1,
		help='When subsampling, first pick rows whose V gene sequences'
			'have at least N exact copies in the input. Default: %(default)s')
	arg('table', help='Table with parsed IgBLAST results')


Groupinfo = namedtuple('Groupinfo', 'count unique_J unique_CDR3')

SiblingInfo = namedtuple('SiblingInfo', 'sequence requested name group')


class SiblingMerger(Merger):
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
		return SiblingInfo(seq, requested, name, group)


class Discoverer:
	"""
	Discover candidates for novel V genes.
	"""
	def __init__(self, database, windows, left, right, cluster, table_output,
			consensus_threshold, v_error_rate, downsample,
			cluster_subsample_size, approx_columns, max_n_bases, exact_copies,
			seed):
		self.database = database
		self.windows = windows
		self.left = left
		self.right = right
		self.cluster = cluster
		self.table_output = table_output
		self.consensus_threshold = consensus_threshold
		self.v_error_rate = v_error_rate
		self.downsample = downsample
		self.cluster_subsample_size = cluster_subsample_size
		self.approx_columns = approx_columns
		self.max_n_bases = max_n_bases
		self.exact_copies = exact_copies
		self.seed = seed

	def _sibling_sequence(self, gene, group):
		"""
		For a given group, compute a consensus sequence over the V gene sequences
		in that group.

		If the found sibling is slightly longer or shorter than the version in
		the database, adjust it so it corresponds to the database version exactly.
		"""
		sequence = iterative_consensus(list(group.V_nt), program='muscle-medium',
			threshold=self.consensus_threshold/100,
			maximum_subsample_size=self.downsample)
		if gene in self.database:
			database_sequence = self.database[gene]
			if sequence.startswith(database_sequence) or database_sequence.startswith(sequence):
				return database_sequence
		return sequence

	@staticmethod
	def _guess_chain(group):
		"""
		Return a guess for the chain type of a given group
		"""
		return Counter(group.chain).most_common()[0][0]

	@staticmethod
	def _guess_cdr3_start(group):
		"""
		Return a guess for the CDR3 start within sequences in the given group
		"""
		return Counter(group.V_CDR3_start).most_common()[0][0]

	def _collect_siblings(self, gene, group):
		"""
		gene -- gene name
		group -- pandas.DataFrame of sequences assigned to that gene
		"""
		group = group.copy()
		# the original reference sequence for all the IgBLAST assignments in this group
		database_sequence = self.database.get(gene, None)
		database_sequence_found = False
		for left, right in self.windows:
			left, right = float(left), float(right)
			group_in_window = group[(left <= group.V_SHM) & (group.V_SHM < right)]
			if len(group_in_window) < MINGROUPSIZE:
				continue
			sibling = self._sibling_sequence(gene, group_in_window)
			if left == int(left):
				left = int(left)
			if right == int(right):
				right = int(right)
			requested = (left, right) == (self.left, self.right)
			name = '{}-{}'.format(left, right)
			if sibling == database_sequence:
				database_sequence_found = True
			yield SiblingInfo(sibling, requested, name, group_in_window)

		if self.cluster:
			if self.exact_copies > 1:
				# Preferentially pick those sequences (for subsampling) that have
				# multiple exact copies, then fill up with the others
				exact_group = group[group.copies >= self.exact_copies]
				indices = downsampled(list(exact_group.index), self.cluster_subsample_size)
				if len(indices) < self.cluster_subsample_size:
					not_exact_group = group[group.copies < self.exact_copies]
					indices.extend(downsampled(list(not_exact_group.index),
						self.cluster_subsample_size - len(indices)))
			else:
				indices = downsampled(list(group.index), self.cluster_subsample_size)
			# Ignore CDR3 part of the V sequence for clustering
			sequences_no_cdr3 = list(group.V_no_CDR3.loc[indices])
			df, linkage, clusters = cluster_sequences(sequences_no_cdr3, MINGROUPSIZE)
			logger.info('Clustering %d sequences (downsampled to %d) assigned to %r gave %d cluster(s)',
				len(group), len(indices), gene, len(set(clusters)))
			cluster_indices = [[] for _ in range(max(clusters) + 1)]
			for i, cluster_id in enumerate(clusters):
				cluster_indices[cluster_id].append(indices[i])

			cl = 1
			for ind in cluster_indices:
				group_in_window = group.loc[ind]
				if len(group_in_window) < MINGROUPSIZE:
					logger.info('Skipping a cluster because it is too small', ind)
					continue
				sibling = self._sibling_sequence(gene, group_in_window)
				name = 'cl{}'.format(cl)
				if sibling == database_sequence:
					database_sequence_found = True
				yield SiblingInfo(sibling, False, name, group_in_window)
				cl += 1

		if database_sequence:
			# If this is a database sequence and there are some exact occurrences of it,
			# add it to the list of candidates even if it has not been found as a
			# cluster.
			group_in_window = group[group.V_errors == 0]
			if len(group_in_window) >= MINGROUPSIZE:
				if not database_sequence_found:
					logger.info('Database sequence seems to be %r expressed, but is missing from '
						'candidates. Re-adding it.', gene)
				yield SiblingInfo(database_sequence, False, 'db', group_in_window)

	def set_random_seed(self, name):
		"""Set random seed depending on gene name and seed given to constructor"""
		h = hashlib.md5(name.encode()).digest()[:4]
		n = int.from_bytes(h, byteorder='big')
		random.seed(n + self.seed)

	def __call__(self, args):
		"""
		Discover new V genes. args is a tuple (gene, group)
		gene -- name of the gene
		group -- a pandas DataFrame with the group corresponding to the gene
		"""
		gene, group = args
		self.set_random_seed(gene)
		siblings = SiblingMerger()
		for sibling in self._collect_siblings(gene, group):
			siblings.add(sibling)

		candidates = []
		for sibling_info in siblings:
			sibling = sibling_info.sequence
			n_bases = sibling.count('N')
			if n_bases > self.max_n_bases:
				logger.debug('Sibling %s has too many N bases', sibling_info.name)
				continue

			sibling_no_cdr3 = sibling[:self._guess_cdr3_start(group)]
			group_exact_V = group[group.V_no_CDR3 == sibling_no_cdr3]
			if self.approx_columns:
				group['consensus_diff'] = [edit_distance(v_no_cdr3, sibling_no_cdr3)
					for v_no_cdr3 in group.V_no_CDR3]
				group_approximate_V = group[group.consensus_diff <= len(sibling_no_cdr3) * self.v_error_rate]
			del sibling_no_cdr3

			groups = (
				('window', sibling_info.group),
				('exact', group_exact_V))
			if self.approx_columns:
				groups += (('approx', group_approximate_V), )
			info = dict()
			for key, g in groups:
				unique_J = len(set(s for s in g.J_gene if s))
				unique_CDR3 = len(set(s for s in g.CDR3_nt if s))
				count = len(g.index)
				info[key] = Groupinfo(count=count, unique_J=unique_J, unique_CDR3=unique_CDR3)
			if gene in self.database:
				database_diff = edit_distance(sibling, self.database[gene])
			else:
				database_diff = None

			# Build the Candidate
			# TODO use UniqueNamer here
			if database_diff is not None and database_diff != 0:
				sequence_id = unique_name(gene, sibling)
			else:
				sequence_id = gene

			if self.approx_columns:
				assert info['exact'].count <= info['approx'].count
				approx = info['approx'].count
				Js_approx = info['approx'].unique_J
				CDR3s_approx = info['approx'].unique_CDR3
			else:
				approx = None
				Js_approx = None
				CDR3s_approx = None

			chain = self._guess_chain(sibling_info.group)
			cdr3_start = self._guess_cdr3_start(sibling_info.group)
			try:
				ratio = info['exact'].count / info['exact'].unique_CDR3
			except ZeroDivisionError:
				ratio = 0
			candidate = Candidate(
				name=sequence_id,
				source=gene,
				chain=chain,
				cluster=sibling_info.name if len(sibling_info.group) < len(group) else 'all',
				cluster_size=info['window'].count,
				Js=info['window'].unique_J,
				CDR3s=info['window'].unique_CDR3,
				exact=info['exact'].count,
				Js_exact=info['exact'].unique_J,
				CDR3s_exact=info['exact'].unique_CDR3,
				CDR3_exact_ratio='{:.2f}'.format(ratio),
				approx=approx,
				Js_approx=Js_approx,
				CDR3s_approx=CDR3s_approx,
				N_bases=n_bases,
				database_diff=database_diff,
				has_stop=int(has_stop(sibling)),
				looks_like_V=int(looks_like_V_gene(sibling, chain)),
				CDR3_start=cdr3_start,
				consensus=sibling,
			)
			candidates.append(candidate)

			# If a window was requested via --left/--right, write the 'approx'
			# subset to a separate file.
			if self.table_output and self.approx_columns and any(si.requested for si in sibling_info) and len(group_approximate_V) > 0:
				if not os.path.exists(self.table_output):
					os.mkdir(self.table_output)
				path = os.path.join(self.table_output, gene + '.tab')
				group_approximate_V.sort('consensus_diff').to_csv(path, sep='\t')
				logger.info('Wrote %s for window %s-%s', path, self.left, self.right)
		return candidates


Candidate = namedtuple('Candidate', [
	'name',
	'source',
	'chain',
	'cluster',
	'cluster_size',
	'Js',
	'CDR3s',
	'exact',
	'Js_exact',
	'CDR3s_exact',
	'CDR3_exact_ratio',
	'approx',
	'Js_approx',
	'CDR3s_approx',
	'N_bases',
	'database_diff',
	'has_stop',
	'looks_like_V',
	'CDR3_start',
	'consensus',
])


def count_prefixes(sequences):
	"""
	Count how often each sequence occurs in the given list of
	sequences. If one sequence is the prefix of another one,
	they are considered to be 'identical'.

	Return a dictionary that maps sequence to count.

	>>> r = count_prefixes(['A', 'BB', 'CD', 'CDE', 'CDEF'])
	>>> r == {'A': 1, 'BB': 1, 'CD': 3, 'CDE': 3, 'CDEF': 3}
	True
	"""
	sequences = sorted(sequences)
	sequences.append('GUARD')
	prev = 'X'
	start = 0
	count = dict()
	for i, s in enumerate(sequences):
		if not s.startswith(prev):
			# end of a run
			for j in range(start, i):
				count[sequences[j]] = i - start
			start = i
		prev = s
	return count


def main(args):
	v_error_rate = args.error_rate / 100
	assert 0 <= v_error_rate <= 1

	if args.database:
		with SequenceReader(args.database) as sr:
			database = {record.name: record.sequence.upper() for record in sr}
	else:
		database = dict()

	if args.seed:
		seed = args.seed
	else:
		seed = random.randrange(10**6)
		logger.info('Use --seed=%d to reproduce this run', seed)

	table = read_table(args.table, usecols=('name', 'chain', 'V_gene', 'J_gene', 'V_nt', 'CDR3_nt', 'V_CDR3_start',
		'V_SHM', 'J_SHM', 'V_errors', 'J_errors'))
	table['V_no_CDR3'] = [s[:start] if start != 0 else s for s, start in
		zip(table.V_nt, table.V_CDR3_start)]

	logger.info('%s rows read', len(table))
	if not args.ignore_J:
		# Discard rows with any mutation within J at all
		table = table[table.J_SHM == 0][:]
		logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	if args.approx:
		logger.info('Approximate comparisons between V gene sequence and consensus '
			'allow %.1f%% errors.', v_error_rate*100)

	if args.exact_copies > 1:
		multiplicities = count_prefixes(table.V_no_CDR3)
		table['copies'] = table.V_no_CDR3.map(multiplicities)
		logger.info('%s rows contain V sequences with at least %s copies',
			sum(table.copies >= args.exact_copies), args.exact_copies)

	columns = list(Candidate._fields)
	if not args.approx:
		for col in ['approx', 'Js_approx', 'CDR3s_approx']:
			columns.remove(col)
	if not args.max_n_bases:
		columns.remove('N_bases')
	writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter='\t', extrasaction='ignore')
	writer.writeheader()
	genes = set(args.gene)
	if args.window_width:
		windows = [(start, start + args.window_width)
			for start in np.arange(0, 20, args.window_width)]
		logger.info('Using an error rate window of %.1f%% to %.1f%%', args.left, args.right)
		windows.append((args.left, args.right))
	else:
		windows = []

	groups = []
	for gene, group in table.groupby('V_gene'):
		if genes and gene not in genes:
			continue
		if len(group) < MINGROUPSIZE:
			continue
		groups.append((gene, group))

	discoverer = Discoverer(database, windows, args.left, args.right, args.cluster,
		args.table_output, args.consensus_threshold, v_error_rate,
		MAXIMUM_SUBSAMPLE_SIZE, cluster_subsample_size=args.subsample,
		approx_columns=args.approx, max_n_bases=args.max_n_bases, exact_copies=args.exact_copies,
		seed=seed)
	n_consensus = 0

	Pool = SerialPool if args.threads == 1 else multiprocessing.Pool
	with Pool(args.threads) as pool:
		for candidates in pool.imap(discoverer, groups, chunksize=1):
			for candidate in candidates:
				writer.writerow(candidate._asdict())
			sys.stdout.flush()
			n_consensus += len(candidates)
	logger.info('%s consensus sequences for %s gene(s) computed', n_consensus, len(groups))
