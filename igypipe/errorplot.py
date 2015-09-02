"""
For each gene, plot a histogram of error rates.

The histogram shows the error rates for all sequences assigned to that gene.
"""
import logging
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from .table import read_table

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('errorplot', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=errorplot_command)
	subparser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene (default: %(default)s)')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')
	subparser.add_argument('pdf', help='Plot error frequency histograms to this PDF file', default=None)
	return subparser


def plot_error_histogram(group, v_gene_name, bins=np.arange(20.1)):
	"""
	Plot a histogram of error rates for a specific V gene.

	v_gene -- name of the gene
	"""
	exact_matches = group[group.V_SHM == 0]
	exact_unique_CDR3 = len(set(s for s in exact_matches.CDR3_nt if s))
	exact_unique_J = len(set(exact_matches.J_gene))

	fig = Figure(figsize=(297/25.4, 210/25.4))
	ax = fig.gca()
	ax.set_xlabel('Error rate (%)')
	ax.set_ylabel('Frequency')
	fig.suptitle('Gene ' + v_gene_name, fontsize=18)
	ax.set_title('{:,} sequences assigned'.format(len(group)), fontsize=16)

	def draw_text(x, y, text):
		ax.text(x, y, text, transform=ax.transAxes, fontsize=16, color='brown', ha='left', va='top')
	draw_text(0.65, 0.95, '{:,} ({:.1%}) exact matches, using'.format(len(exact_matches), len(exact_matches) / len(group)))
	draw_text(0.7, 0.9, '{} unique CDR3s'.format(exact_unique_CDR3))
	draw_text(0.7, 0.85, '{} unique Js'.format(exact_unique_J))

	_ = ax.hist(list(group.V_SHM), bins=bins)
	return fig


def errorplot_command(args):
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	with PdfPages(args.pdf) as pages:
		for gene, group in table.groupby('V_gene'):
			if len(group) < args.minimum_group_size:
				continue
			fig = plot_error_histogram(group, gene)
			n += 1
			FigureCanvasPdf(fig).print_figure(pages)
	logger.info('%s plots created (rest had too few sequences)', n)
