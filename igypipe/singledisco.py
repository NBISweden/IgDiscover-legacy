"""
Discover new V genes within a single antibody library.
"""
import logging
from collections import Counter, OrderedDict
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from sqt.align import multialign, consensus, edit_distance
from .table import read_table

logger = logging.getLogger(__name__)

MINGROUPSIZE_PLOT = 200
MINGROUPSIZE_CONSENSUS = 10


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--plot', help='Plot error frequency histograms to this file', default=None)
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Use "all" to compute for all genes.')
	subparser.add_argument('--left', '-l', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at least this %%SHM (default: %(default)s)', default=0)
	subparser.add_argument('--right', '-r', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at most this %%SHM (default: %(default)s)', default=100)
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


def sister_sequence(group, min_shm, max_shm, program='muscle-medium'):
	"""

	"""
	sequences = OrderedDict()
	# TODO Perhaps create the dict in such a way that those with the most
	# abundant no. of errors are put in first.
	for _, row in group.iterrows():
		if min_shm <= row.V_SHM <= max_shm:
			sequences[row.name] = row.V_nt
	logger.info('Computing consensus from %s sequences', len(sequences))
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
				if len(group) < MINGROUPSIZE_PLOT:
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
		if len(group) < MINGROUPSIZE_CONSENSUS:
			logger.info('Skipping %s as the number of sequences is too small (%s)', gene, len(group))
			continue
		s = sister_sequence(group, min_shm=args.left, max_shm=args.right)
		print('>{}_sister\n{}'.format(gene, s))
		n_count = s.count('N')
		exact_count = sum(group.V_nt == s)
		distances = [ edit_distance(v_nt, s) for v_nt in group.V_nt ]
		less_than_one = sum(d <= len(s) * 0.01 for d in distances)
		logger.info('Consensus for %s has %s “N” bases', gene, n_count)
		logger.info('There are %s exact occurrences in all input sequences', exact_count)
		logger.info('There are %s (%.1f%%) occurrences with at most 1%% differences',
			  less_than_one, less_than_one / len(distances) * 100)
		n += 1
	if genes and n > 1:
		logger.info('%s consensus sequences computed', n)
