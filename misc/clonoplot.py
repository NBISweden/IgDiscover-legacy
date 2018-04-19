"""
Query two ...
"""
import logging
from itertools import islice
from collections import Counter
from io import StringIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
import pandas as pd
import seaborn as sns

sns.set()#style='white')
logger = logging.getLogger(__name__)
CM = 1 / 2.54


def add_arguments(parser):
	arg = parser.add_argument
	arg('--nt', default=False, action='store_true',
		help='Count CDR3 mismatches on nucleotide level. Default: Compare amino-acids.')
	arg('--limit', default=None, type=int, metavar='N',
		help='Read only the first N clonotypes of each table (for fast testing)')
	arg('--minsize', default=5, type=int, metavar='N',
		help='Require at least N members in each dataset to plot a clonotype')
	arg('pdf', help='PDF output')
	arg('tables', nargs='+', help='clonotype member tables')


def read_dataset(path, cdr3_column='CDR3_aa', limit=None):
	usecols = [
		'V_gene',
		'J_gene',
		cdr3_column,
		'FR1_SHM',
		'CDR1_SHM',
		'FR2_SHM',
		'CDR2_SHM',
		'FR3_SHM',
		'V_SHM',
		'J_SHM',
	]
	column_names = list(pd.read_csv(path, sep='\t', nrows=0).columns)
	assert 'CDR3_length' in column_names
	with open(path) as f:
		data = f.read()
	first = True
	tables = []
	for chunk in islice(data.split('\n\n'), limit):
		table = pd.read_csv(StringIO(chunk), sep='\t', usecols=usecols, header=0 if first else None,
			names=column_names)
		first = False
		assert len(set(table[cdr3_column])) == 1
		tables.append(table.set_index(['V_gene', 'J_gene', cdr3_column]).sort_index())
	return pd.concat(tables)


def main(args):
	cdr3_column = 'CDR3_nt' if args.nt else 'CDR3_aa'
	n_datasets = len(args.tables)
	datasets = (read_dataset(path, limit=args.limit) for path in args.tables)
	df = pd.concat(datasets, keys=range(n_datasets), names=['dataset_id'])
	logger.info('Read %s tables', n_datasets)
	df.rename(columns=lambda x: x[:-4], inplace=True)
	cols = ['V_gene', 'J_gene', cdr3_column]
	n = 0
	with PdfPages(args.pdf) as pages:
		for (v_gene, j_gene, cdr3), group in df.groupby(level=cols):
			group = group.reset_index(level=cols, drop=True)
			skip = False
			counter = Counter(group.index)
			for dataset_id in range(n_datasets):
				if counter[dataset_id] < args.minsize:
					skip = True
					break
			if skip:
				continue
			table = group.stack()
			table.index.set_names('region', level=1, inplace=True)
			table.name = 'SHM'
			g = sns.factorplot(data=table.reset_index(), x='region', y='SHM', hue='dataset_id',
				kind='violin', size=5, aspect=2)
			dscounts = ' vs '.join(str(counter[i]) for i in range(n_datasets))
			g.fig.suptitle('V: {} – J: {} – CDR3: {} ({})'.format(v_gene, j_gene, cdr3, dscounts))
			g.set_axis_labels('Region')
			g.set_ylabels('%SHM')
			# g.despine()
			FigureCanvasPdf(g.fig).print_figure(pages, bbox_inches='tight')
			n += 1
			logger.info('Plotted %s clonotypes', n)


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
