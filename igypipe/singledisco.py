"""
Discover new V genes within a single antibody library.
"""
import logging
from collections import Counter, OrderedDict
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from sqt.align import multialign, consensus
from .table import read_table

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--plot', help='Plot error frequency histograms to this file', default=None)
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Use "all" to compute for all genes.')
	subparser.add_argument('--threshold', type=float, help='Compute consensus from those sequences that have at least this %%SHM (default: %(default)s)', default=0)

	#subparser.add_argument('--minimum-frequency', '-n', type=int, metavar='N',
		#default=None,
		#help='Minimum number of datasets in which sequence must occur (default is no. of files divided by two)')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser



def plot_shms(group, v_gene, bins=np.arange(20)):
	"""
	Plot error frequency distribution for a specific V gene.

	v_gene -- name of the gene
	"""
	shms = list(group.V_SHM)
	#mean = np.mean(shms)
	z = shms.count(0)
	fig = Figure(figsize=(297/25.4, 210/25.4))
	ax = fig.gca()
	ax.set_xlabel('%SHM')
	ax.set_ylabel('Frequency')
	ax.set_title('Gene ' + v_gene, fontsize=18)
	ax.text(0.95, 0.95, 'zero differences: {} times'.format(z), transform=ax.transAxes, fontsize=15, ha='right', va='top')

	#ax.axvline(mean, color='darkred')
	_ = ax.hist(shms, bins=bins)
	return fig


def sister_sequence(group, shm_threshold, program='muscle-medium'):
	"""

	"""
	sequences = OrderedDict()
	# TODO Perhaps create the dict in such a way that those with the most
	# abundant no. of errors are put in first.
	for _, row in group.iterrows():
		if row.V_SHM >= shm_threshold:
			sequences[row.name] = row.V_nt
	aligned = multialign(sequences, program=program)
	return consensus(aligned, threshold=0.6)


def discover_command(args):
	table = read_table(args.table)
	#table = table.loc[:,['V_gene', 'V%SHM', 'V_nt', 'name']]
	#t = table.loc[:,('name', 'V_gene', 'V_nt', 'V%SHM')]

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	if args.plot:
		n = 0
		with PdfPages(args.plot) as pages:
			for gene, group in table.groupby('V_gene'):
				if len(group) < 200:
					continue
				fig = plot_shms(group, gene)
				n += 1
				FigureCanvasPdf(fig).print_figure(pages)
		logger.info('%s plots created (rest had too few sequences)', n)

	genes = set(args.gene)
	n = 0
	for gene, group in table.groupby('V_gene'):
		if not ('all' in genes or gene in genes):
			continue
		if len(group) < 10:
			continue
		s = sister_sequence(group, args.threshold)
		print('>{}_sister\n{}'.format(gene, s))
		n += 1
	logger.info('%s consensus sequences computed', n)
