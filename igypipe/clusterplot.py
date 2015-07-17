"""
For each V gene, plot a clustermap ...
"""
import sys
import logging
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sqt.align import edit_distance
from .table import read_table
from .utils import downsampled

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('clusterplot', help=__doc__)
	subparser.set_defaults(func=command)
	subparser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene (default: %(default)s)')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	subparser.add_argument('pdf', help='Plot clustermaps to this PDF file', default=None)
	return subparser


def distances(sequences, band=0.2):
	"""
	Compute all pairwise edit distances and return a square matrix.

	Entry [i,j] in the matrix is the edit distance between sequences[i]
	and sequences[j].
	"""
	m = np.zeros((len(sequences), len(sequences)), dtype=float)
	maxdiff = max(int(len(s) * band) for s in sequences)
	for i, s in enumerate(sequences):
		for j, t in enumerate(sequences):
			if i < j:
				m[j,i] = m[i,j] = min(maxdiff+1, edit_distance(s, t, maxdiff=maxdiff))
	return m


def find_best(tree):
	if tree.count == 1:
		return None
	if tree.left.count <= 2:
		assert tree.right is not None
		return find_best(tree.right)
	if tree.right.count <= 2:
		assert tree.left is not None
		return find_best(tree.left)
	return tree


def collect_ids(tree):
	if tree.count == 1:
		return [tree.id]
	return collect_ids(tree.left) + collect_ids(tree.right)


def plot_clustermap(group, v_gene):
	"""
	Plot a clustermap
	for a specific V gene.

	v_gene -- name of the gene
	"""
	sequences = list(group.V_nt)
	sequences = downsampled(sequences, 300)
	m = distances(sequences)
	d = pd.DataFrame(m)
	dists = distance.squareform(m)
	linkage = hierarchy.linkage(dists, method='average')

	palette = sns.color_palette(['red', 'blue', 'black'])

	##clusters = hierarchy.fcluster(linkage, 0.5*(linkage[-2,2] + linkage[-1,2]), criterion='distance')
	##print('clusters:', clusters)
	n = find_best(hierarchy.to_tree(linkage))
	if n is not None:
		left_ids, right_ids = collect_ids(n.left), collect_ids(n.right)
		#print(n.id, n.left.count, n.right.count)
	else:
		left_ids, right_ids = [], []
	row_colors = [palette[0] if i in left_ids else palette[1] if i in right_ids else palette[2] for i in range(len(d))]
	#col_colors = [palette[1-c] for c in clusters]
	cm = sns.clustermap(d,
			row_colors=row_colors, row_linkage=linkage,
			col_linkage=linkage,
			linewidths=0, linecolor='none', figsize=(210/25.4, 210/25.4), cmap='Blues',
			xticklabels=False, yticklabels=False)
	return cm.fig


def command(args):
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	try:
		with PdfPages(args.pdf) as pages:
			for gene, group in table.groupby('V_gene'):
				if len(group) < args.minimum_group_size:
					continue
				logger.info('Working on %s with %s sequences', gene, len(group))
				fig = plot_clustermap(group, gene)
				fig.suptitle(gene)
				FigureCanvasPdf(fig).print_figure(pages)
				n += 1
				if n == 4:
					break
	except KeyboardInterrupt:
		logger.warn('Cancelled')
		sys.exit(130)
	logger.info('%s plots created (rest had too few sequences)', n)
