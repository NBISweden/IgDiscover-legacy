"""
Plot histograms of differences to reference V gene

For each gene, a histogram is plotted that shows how often a sequence was
assigned to that gene at a certain percentage difference.
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
	parser.add_argument('--ignore-J', action='store_true', default=False,
		help='Include also rows without J assignment or J%%SHM>0.')
	parser.add_argument('table', help='Table with parsed IgBLAST results')
	parser.add_argument('pdf', help='Plot error frequency histograms to this PDF file', default=None)


def plot_difference_histogram(group, gene_name, bins=np.arange(20.1)):
	"""
	Plot a histogram of percentage differences for a specific gene.
	"""
	exact_matches = group[group.V_SHM == 0]
	CDR3s_exact = len(set(s for s in exact_matches.CDR3_nt if s))
	Js_exact = len(set(exact_matches.J_gene))

	fig = Figure(figsize=(100/25.4, 60/25.4))
	ax = fig.gca()
	ax.set_xlabel('Percentage difference')
	ax.set_ylabel('Frequency')
	fig.suptitle('Gene ' + gene_name, y=1.08, fontsize=16)
	ax.set_title('{:,} sequences assigned'.format(len(group)))

	ax.text(0.25, 0.95,
		'{:,} ({:.1%}) exact matches\n  {} unique CDR3\n  {} unique J'.format(
			len(exact_matches), len(exact_matches) / len(group),
			CDR3s_exact, Js_exact),
		transform=ax.transAxes, fontsize=10,
		bbox=dict(boxstyle='round', facecolor='white', alpha=0.5),
		horizontalalignment='left', verticalalignment='top')
	_ = ax.hist(list(group.V_SHM), bins=bins)
	return fig


def main(args):
	table = read_table(args.table)

	# Discard rows with any mutation within J at all
	logger.info('%s rows read', len(table))
	if not args.ignore_J:
		# Discard rows with any mutation within J at all
		table = table[table.J_SHM == 0][:]
		logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	n = 0
	too_few = 0
	with PdfPages(args.pdf) as pages:
		for gene, group in table.groupby('V_gene'):
			if len(group) < args.minimum_group_size:
				too_few += 1
				continue
			fig = plot_difference_histogram(group, gene)
			n += 1
			FigureCanvasPdf(fig).print_figure(pages, bbox_inches='tight')
	logger.info('%s plots created (%s skipped because of too few sequences)', n, too_few)
