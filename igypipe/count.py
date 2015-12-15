"""
Count and plot V, D, J gene usage.
"""
import sys
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
	table = read_table(args.table, log=True)

	columns = ('gene', 'count', 'unique_CDR3')
	if args.gene != 'V':
		columns += ('unique_V',)
	if args.gene != 'J':
		columns += ('unique_J',)
	column_name = '{}_gene'.format(args.gene)
	if len(table) > 0:
		rows = []
		groups = table.groupby(column_name)
		for gene, group in groups:
			count = len(group)
			unique_CDR3 = len(set(s for s in group.CDR3_nt if s))
			row = dict(gene=gene, count=count, unique_CDR3=unique_CDR3)
			if args.gene != 'V':
				row['unique_V'] = len(set(s for s in group.V_gene if s))
			if args.gene != 'J':
				row['unique_J'] = len(set(s for s in group.J_gene if s))
			rows.append(row)
		counts = pd.DataFrame(rows, columns=columns).set_index('gene')
	else:
		# Work around a pandas bug in reindex when the table is empty
		counts = pd.DataFrame(columns=columns, dtype=int).set_index('gene')

	# Make sure that always all gene names are listed even if no sequences
	# were assigned.
	if gene_names:
		counts = counts.reindex(gene_names, fill_value=0)

	counts.to_csv(sys.stdout, sep='\t')
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
	ax.set_xlim((-0.25, len(counts)))
	ax.bar(np.arange(len(counts)), counts['count'])
	fig.set_tight_layout(True)
	fig.savefig(args.plot)
	logger.info("Wrote %s", args.plot)
