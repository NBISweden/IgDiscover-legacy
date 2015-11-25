"""
Count and plot V, D, J gene usage.
"""
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from sqt import FastaReader

from .table import read_table
from .utils import natural_sort_key

logger = logging.getLogger(__name__)

def add_arguments(parser):
	arg = parser.add_argument
	arg('--gene', default='V', choices=('V', 'D', 'J'),
		help='Which type gene to count. Default: V')
	arg('--database', metavar='FASTA',
		help='FASTA file with reference gene sequences. The names are used to '
		'ensure all names appear in the plot')
	arg('table', help='Table with parsed and filtered IgBLAST results.')
	arg('plot', nargs='?', help='Plot file (png or pdf).')


def main(args):
	if args.database:
		with FastaReader(args.database) as fr:
			gene_names = [ record.name for record in fr ]
		gene_names.sort(key=natural_sort_key)
	else:
		gene_names = None
	d = read_table(args.table, log=True)

	column_name = '{}_gene'.format(args.gene)
	# Work around a pandas bug in reindex when the table is empty
	if len(d) > 0:
		counts = d.groupby(column_name).size()
	else:
		counts = pd.Series([], dtype=int)

	# Make sure that always all gene names are listed, even if expression is 0.
	if gene_names:
		counts = counts.reindex(gene_names, fill_value=0)

	print(counts.to_csv(None, sep='\t'))
	logger.info("Wrote CSV")

	if not args.plot:
		return

	fig = plt.figure(figsize=((50 + len(counts) * 5) / 25.4, 210/25.4))
	matplotlib.rcParams.update({'font.size': 14})
	ax = fig.gca()
	ax.set_title('{} gene usage'.format(args.gene))
	ax.set_xlabel('{} gene'.format(args.gene))
	ax.set_ylabel('Count')
	ax.set_xticks(np.arange(len(counts)) + 0.5)
	ax.set_xticklabels(counts.index, rotation='vertical')
	ax.grid(axis='x')
	ax.set_xlim((-0.5, None))
	ax.bar(np.arange(len(counts)), counts)
	fig.tight_layout()
	fig.savefig(args.plot)
	logger.info("Wrote %s", args.plot)
