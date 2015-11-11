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


def add_arguments(parser):
	parser.add_argument('--minimum-group-size', '-m', metavar='N', default=200,
		help='Do not plot if there are less than N sequences for a gene (default: %(default)s)')
	parser.add_argument('table', help='Table with parsed IgBLAST results')
	parser.add_argument('pdf', help='Plot error frequency histograms to this PDF file', default=None)


def plot_error_histogram(group, v_gene_name, bins=np.arange(20.1)):
	"""
	Plot a histogram of error rates for a specific V gene.

	v_gene -- name of the gene
	"""
	exact_matches = group[group.V_SHM == 0]
	exact_unique_CDR3 = len(set(s for s in exact_matches.CDR3_nt if s))
	exact_unique_J = len(set(exact_matches.J_gene))

	fig = Figure(figsize=(297/2/25.4, 210/2/25.4))
	ax = fig.gca()
	ax.set_xlabel('Error rate (%)')
	ax.set_ylabel('Frequency')
	fig.suptitle('Gene ' + v_gene_name, y=1.02, fontsize=16)
	ax.set_title('{:,} sequences assigned'.format(len(group)))

	ax.text(0.3, 0.95,
		'{:,} ({:.1%}) exact matches, using\n  {} unique CDR3\n  {} unique J'.format(
			len(exact_matches), len(exact_matches) / len(group),
			exact_unique_CDR3, exact_unique_J),
		transform=ax.transAxes, fontsize=14,
		bbox=dict(boxstyle='round', facecolor='white', alpha=0.5),
		horizontalalignment='left', verticalalignment='top')
	_ = ax.hist(list(group.V_SHM), bins=bins)
	return fig


def main(args):
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read', len(table))
	table = table[table.J_SHM == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	with PdfPages(args.pdf) as pages:
		for gene, group in table.groupby('V_gene'):
			if len(group) < args.minimum_group_size:
				continue
			fig = plot_error_histogram(group, gene)
			n += 1
			FigureCanvasPdf(fig).print_figure(pages, bbox_inches='tight')
	logger.info('%s plots created (rest had too few sequences)', n)
