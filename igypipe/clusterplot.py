"""
For each V gene, plot a clustermap of the sequences assigned to it.
"""
import os.path
import logging
import seaborn as sns
from .table import read_table
from .utils import downsampled
from .cluster import cluster_sequences

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene. Default: %(default)s')
	arg('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Default: Compute for all genes.')
	arg('--size', metavar='N', type=int, default=300,
		help='Show at most N sequences (with a matrix of size N x N). Default: %(default)s')
	arg('table', help='Table with parsed and filtered IgBLAST results')
	arg('directory', help='Save clustermaps as PNG into this directory', default=None)


def plot_clustermap(group, gene, plotpath, size=300):
	"""
	Plot a clustermap for a specific V gene.

	size -- Downsample to this many sequences
	gene -- gene name (only used to plot the title)

	Return the number of clusters.
	"""
	sequences = list(group.V_nt)
	sequences = downsampled(sequences, size)
	df, linkage, clusters = cluster_sequences(sequences)

	palette = sns.color_palette(['black']) + sns.color_palette('Set1', n_colors=max(clusters), desat=.8)
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

	return len(set(clusters))


def main(args):
	if not os.path.exists(args.directory):
		os.mkdir(args.directory)
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	genes = frozenset(args.gene)
	n = 0
	too_few = 0
	for gene, group in table.groupby('V_gene'):
		if genes and gene not in genes:
			continue
		if len(group) < args.minimum_group_size:
			too_few += 1
			continue
		n_clusters = plot_clustermap(group, gene, os.path.join(args.directory, gene + '.png'), size=args.size)
		n += 1
		logger.info('Plotted %r with %d clusters', gene, n_clusters)
		#for i, cons in enumerate(consensus_sequences):
			#print('>{}_cluster{}\n{}'.format(gene, ('red', 'blue')[i], cons))
			#print('number of Ns:', cons.count('N'))
		#if len(consensus_sequences) >= 2:
			#print('difference between consensuses:', edit_distance(*consensus_sequences[:2]))


	logger.info('%s plots created (%s skipped because of too few sequences)', n, too_few)
