"""
Discover new V genes within a single antibody library.
"""
import logging
from collections import Counter
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
import seaborn as sns
import numpy as np
from .table import read_table

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('singledisco', help=__doc__)
	subparser.set_defaults(func=discover_command)
	subparser.add_argument('--plot', help='Plot error frequency histograms to this file', default=None)
	#subparser.add_argument('--minimum-frequency', '-n', type=int, metavar='N',
		#default=None,
		#help='Minimum number of datasets in which sequence must occur (default is no. of files divided by two)')
	subparser.add_argument('table', help='Table with parsed IgBLAST results')  # nargs='+'
	return subparser



def plot_shms(table, v_gene, bins=np.arange(20)):
	"""
	Plot error frequency distribution for a specific V gene.

	v_gene -- name of the gene
	"""
	shms = list(table[table.V_gene == v_gene]['V%SHM'])
	if len(shms) < 200:
		return None
	#mean = np.mean(shms)
	z = shms.count(0)
	fig = Figure(figsize=(297/25.4, 210/25.4))
	ax = fig.gca()
	ax.set_xlabel('%SHM')
	ax.set_ylabel('Frequency')
	ax.set_title('Gene ' + v_gene, fontsize=18)
	ax.text(0.95, 0.95, 'zero differences: {} times'.format(z), transform=ax.transAxes, fontsize=15, ha='right', va='top')

	#ax.axvline(mean, color='darkred')
	_ = ax.hist(list(shms), bins=bins)
	return fig


def discover_command(args):
	table = read_table(args.table)
	#table = table.loc[:,['V_gene', 'V%SHM', 'V_nt', 'name']]
	#t = table.loc[:,('name', 'V_gene', 'V_nt', 'V%SHM')]

	# Discard rows with any mutation within J at all
	logger.info('%s rows read (filtered)', len(table))
	table = table[table['J%SHM'] == 0][:]
	logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

	if args.plot:
		n = 0
		with PdfPages(args.plot) as pages:
			for gene in set(table.V_gene):
				fig = plot_shms(table, gene)
				if fig is None:
					continue
				n += 1
				FigureCanvasPdf(fig).print_figure(pages)
		logger.info('%s plots created (rest had too few sequences)', n)
	return


	# Count V sequence occurrences
	counter = Counter()
	for table in tables:
		counter.update(set(table.V_nt))

	# Find most frequent occurrences and print result
	print('Frequency', 'Gene', '%SHM', 'Sequence', sep='\t')
	for sequence, frequency in counter.most_common():
		if frequency < minimum_frequency:
			break
		names = []
		gene = None
		for table in tables:
			matching_rows = table[table.V_nt == sequence]
			if matching_rows.empty:
				continue
			names.extend(matching_rows.name)
			if gene is None:
				row = matching_rows.iloc[0]
				gene = row['V_gene']
				shm = row['V%SHM']
		print(frequency, gene, shm, sequence, *names, sep='\t')
