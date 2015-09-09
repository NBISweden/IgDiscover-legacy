"""
For each V gene, plot a clustermap of the sequences assigned to it.
"""
import os.path
import logging
import seaborn as sns
from .table import read_filtered_table
from .utils import downsampled
from .cluster import cluster_sequences

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('clusterplot', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=command)
	subparser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene. Default: %(default)s')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	subparser.add_argument('directory', help='Save clustermaps as PNG into this directory', default=None)
	return subparser


def plot_clustermap(group, gene, plotpath):
	"""
	Plot a clustermap for a specific V gene.

	gene -- gene name (only used to plot the title)
	"""
	sequences = list(group.V_nt)
	sequences = downsampled(sequences, 300)
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


def command(args):
	if not os.path.exists(args.directory):
		os.mkdir(args.directory)
	table = read_filtered_table(args.table)

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
