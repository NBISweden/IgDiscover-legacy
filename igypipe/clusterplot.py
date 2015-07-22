"""
For each V gene, plot a clustermap of the sequences assigned to it.
"""
import sys
import os.path
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sqt.utils import available_cpu_count
from sqt.align import edit_distance
from .table import read_table
from .utils import downsampled, distances

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('clusterplot', help=__doc__)
	subparser.set_defaults(func=command)
	subparser.add_argument('--threads', '-j', type=int, default=2,
		help='Number of threads. Default: %(default)s')
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


def collect_ids(root):
	"""
	Return a list of ids of all leaves of the given tree
	"""
	if root.is_leaf():
		return [root.id]
	return collect_ids(root.left) + collect_ids(root.right)


def plot_clustermap(group, gene, pdfpath):
	"""
	Plot a clustermap for a specific V gene.

	gene -- gene name (only used to plot the title)
	"""
	sequences = list(group.V_nt)
	sequences = downsampled(sequences, 300)
	m = distances(sequences)
	d = pd.DataFrame(m)
	dists = distance.squareform(m)
	linkage = hierarchy.linkage(dists, method='average')

	palette = sns.color_palette(['red', 'blue', 'black'])

	n = find_separating_node(hierarchy.to_tree(linkage))
	if n is not None:
		left_ids, right_ids = collect_ids(n.left), collect_ids(n.right)
	else:
		left_ids, right_ids = [], []
	row_colors = [palette[0] if i in left_ids else palette[1] if i in right_ids else palette[2] for i in range(len(d))]
	cm = sns.clustermap(d,
			row_linkage=linkage,
			col_linkage=linkage,
			row_colors=row_colors,
			linewidths=0, linecolor='none', figsize=(210/25.4, 210/25.4), cmap='Blues',
			xticklabels=False, yticklabels=False)
	cm.fig.suptitle(gene)
	cm.savefig(pdfpath, dpi=200)
	return gene


def command(args):
	if not os.path.exists(args.directory):
		os.mkdir(args.directory)

	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	with ProcessPoolExecutor(max_workers=args.threads) as executor:
		futures = []
		for gene, group in table.groupby('V_gene'):
			if len(group) < args.minimum_group_size:
				continue
			futures.append(executor.submit(plot_clustermap, group, gene, os.path.join(args.directory, gene + '.png')))
			n += 1
		for future in as_completed(futures):
			gene = future.result()
			logger.info('Plotted %r', gene)

	logger.info('%s plots created (rest had too few sequences)', n)
