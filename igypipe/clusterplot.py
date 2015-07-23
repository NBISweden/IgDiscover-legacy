"""
For each V gene, plot a clustermap of the sequences assigned to it.
"""
import sys
import os.path
import logging
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sqt.align import edit_distance
from .table import read_table
from .utils import downsampled, distances, iterative_consensus

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('clusterplot', help=__doc__)
	subparser.set_defaults(func=command)
	subparser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene. Default: %(default)s')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	subparser.add_argument('directory', help='Save clustermaps as PNG into this directory', default=None)
	return subparser


def inner_nodes(root):
	"""
	Return a list of all inner nodes of the tree
	"""
	if root.is_leaf():
		return []
	return inner_nodes(root.left) + [root] + inner_nodes(root.right)


def find_separating_node_old(root):
	if root.count == 1:
		return None
	if root.left.count <= 2:
		return find_best(root.right)
	if root.right.count <= 2:
		return find_best(root.left)
	return root


def find_separating_node(root):
	"""
	Return an inner node of the tree such that the left and right subtrees
	of that node describe different clusters. A simple heuristic is used
	that inspects all inner nodes and skips those where one of the two joined
	subtrees is a singleton. Then for each node the size of the smaller subtree
	is computed. Finally, the node where that value is largest is picked.
	"""
	best = None
	best_count = 0
	for node in inner_nodes(root):
		if node.left.count == 1 or node.right.count == 1:
			continue
		smallest_count = min(node.left.count, node.right.count)
		if smallest_count > best_count:
			best = node
			best_count = smallest_count
	return best


def find_clusters(root):
	candidates = []
	for node in inner_nodes(root):
		if node.left.count == 1 or node.right.count == 1:
			continue
		smallest_count = min(node.left.count, node.right.count)

		if smallest_count > best_count:
			best = node
			best_count = smallest_count
	return best


def collect_ids(root):
	"""
	Return a list of ids of all leaves of the given tree
	"""
	if root.is_leaf():
		return [root.id]
	return collect_ids(root.left) + collect_ids(root.right)


def cluster_sequences_old(sequences):
	"""
	Cluster the given sequences into groups of similar sequences.

	Return a triple that contains a pandas.DataFrame with the edit distances,
	the linkage result, and a list that maps sequence ids to their cluster id.
	If an entry is zero in that list, it means that the sequence is not part of
	a cluster.
	"""
	matrix = distances(sequences)
	linkage = hierarchy.linkage(distance.squareform(matrix), method='average')
	sepnode = find_separating_node(hierarchy.to_tree(linkage))
	id_lists = []
	if sepnode is not None:
		for n in (sepnode.left, sepnode.right):
			id_lists.append(collect_ids(n))
	clusters = [0] * len(sequences)
	for i, ids in enumerate(id_lists, start=1):
		for id in ids:
			clusters[id] = i

	return pd.DataFrame(matrix), linkage, clusters


def cluster_sequences(sequences):
	"""
	Cluster the given sequences into groups of similar sequences.

	Return a triple that contains a pandas.DataFrame with the edit distances,
	the linkage result, and a list that maps sequence ids to their cluster id.
	If an entry is zero in that list, it means that the sequence is not part of
	a cluster.
	"""
	matrix = distances(sequences)
	linkage = hierarchy.linkage(distance.squareform(matrix), method='average')
	inner = inner_nodes(hierarchy.to_tree(linkage))
	prev = 100000
	clusters = [0] * len(sequences)
	cl = 1
	for n in inner:
		if prev/n.dist < 0.8 and n.left.count >= 5 and n.right.count >= 5:
			for id in collect_ids(n.left):
				clusters[id] = cl
			cl += 1
		prev = n.dist

	return pd.DataFrame(matrix), linkage, clusters


def cluster_consensus(sequences, clusters, minsize=5):
	"""
	Compute a consensus for each cluster that has at least minsize members.

	Return a list of consensus sequences.
	"""
	n_clusters = max(clusters) + 1
	cluster_sequences = [ [] for _ in range(n_clusters) ]
	for i, cluster_id in enumerate(clusters):
		if cluster_id == 0:
			# sequence not assigned to a cluster
			continue
		cluster_sequences[cluster_id-1].append(sequences[i])
	cons = []
	for seqs in cluster_sequences:
		if len(seqs) >= minsize:
			cons.append(iterative_consensus(seqs))
	return cons


def plot_clustermap(group, gene, plotpath):
	"""
	Plot a clustermap for a specific V gene.

	gene -- gene name (only used to plot the title)
	"""
	sequences = list(group.V_nt)
	sequences = downsampled(sequences, 300)
	df, linkage, clusters = cluster_sequences(sequences)

	palette = sns.color_palette(['black']) + sns.color_palette('Set1', n_colors=20, desat=.8)
	row_colors = [ palette[cluster_id] for cluster_id in clusters ]
	cm = sns.clustermap(df,
			row_linkage=linkage,
			col_linkage=linkage,
			row_colors=row_colors,
			linewidths=0, linecolor='none', figsize=(210/25.4, 210/25.4), cmap='Blues',
			xticklabels=False, yticklabels=False)
	cm.fig.suptitle(gene)
	cm.savefig(plotpath, dpi=200)

	# free the memory used by the plot
	import matplotlib.pyplot as plt
	plt.close('all')


def command(args):
	if not os.path.exists(args.directory):
		os.mkdir(args.directory)
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	for gene, group in table.groupby('V_gene'):
		if len(group) < args.minimum_group_size:
			continue
		plot_clustermap(group, gene, os.path.join(args.directory, gene + '.png'))
		n += 1
		logger.info('Plotted %r', gene)
		#for i, cons in enumerate(consensus_sequences):
			#print('>{}_cluster{}\n{}'.format(gene, ('red', 'blue')[i], cons))
			#print('number of Ns:', cons.count('N'))
		#if len(consensus_sequences) >= 2:
			#print('difference between consensuses:', edit_distance(*consensus_sequences[:2]))


	logger.info('%s plots created (rest had too few sequences)', n)
