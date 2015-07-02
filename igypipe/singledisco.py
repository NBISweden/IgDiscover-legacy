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
	subparser.add_argument('--error-rate', metavar='PERCENT', type=float, default=1,
		help='When finding approximate V gene matches, allow PERCENT errors (default: %(default)s)')
	subparser.add_argument('--gene', '-g', action='append', default=[],
		help='Compute consensus for this gene. Can be given multiple times. Use "all" to compute for all genes.')
	subparser.add_argument('--left', '-l', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at least this %%SHM (default: %(default)s)', default=0)
	subparser.add_argument('--right', '-r', type=float, metavar='%SHM',
		help='For consensus, include only sequences that have at most this %%SHM (default: %(default)s)', default=100)
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser


def plot_shms(group, v_gene, bins=np.arange(20.1)):
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
	ax.text(0.95, 0.95, '{} sequences with zero differences'.format(z), transform=ax.transAxes, fontsize=15, ha='right', va='top')
	ax.text(0.95, 0.90, '{} different J genes used'.format(len(set(group.J_gene))), transform=ax.transAxes, fontsize=15, ha='right', va='top')

	#ax.axvline(mean, color='darkred')
	_ = ax.hist(shms, bins=bins)
	return fig


def sister_sequence(group, program='muscle-medium'):
	"""

	"""
	sequences = OrderedDict()
	# TODO Perhaps create the dict in such a way that those with the most
	# abundant no. of errors are put in first.
	for _, row in group.iterrows():
		sequences[row.name] = row.V_nt
	logger.info('Computing consensus from %s sequences', len(sequences))
	aligned = multialign(sequences, program=program)
	return consensus(aligned, threshold=0.6)


def discover_command(args):
	v_error_rate = args.error_rate / 100
	assert 0 <= v_error_rate <= 1
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

		# Create various subsets of the full group
		group_in_shm_range = group[(group.V_SHM >= args.left) & (group.V_SHM <= args.right)]
		s = sister_sequence(group_in_shm_range)
		logger.info('Consensus for %s has %s “N” bases', gene, s.count('N'))
		group_exact_V = group[group.V_nt == s]
		group_approximate_V = group[list(edit_distance(v_nt, s) <= len(s) * v_error_rate for v_nt in group.V_nt)]

		for description, g in (
				('sequences in total', group),
				('sequences are within the given %SHM range', group_in_shm_range),
				('sequences match the consensus exactly', group_exact_V),
				('sequences match the consens approximately', group_approximate_V)):
			logger.info('%s %s (%.1f%%):', len(g), description, len(g) / len(group) * 100)
			logger.info('   %s different J genes', len(set(g.J_gene)))
			logger.info('   %s different CDR3 sequences', len(set(g.CDR3_nt)))

		# Print this last so it doesn’t mess up output too bad in case stdout
		# isn’t redirected anywhere.
		print('>{}_sister\n{}'.format(gene, s))
		n += 1
	if genes and n > 1:
		logger.info('%s consensus sequences computed', n)
