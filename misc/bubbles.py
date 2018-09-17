#!/usr/bin/env python3
"""
Bubble plot
"""
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#sns.set(style='white', font_scale=1.5, rc={"lines.linewidth": 1})
logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('pdf', help='PDF output')
	arg('table', help='Input table')


def main(args):
	df = pd.read_table(args.table, sep=';', index_col=0)
	dfu = df.unstack().reset_index()
	dfu.columns = ['compartment', 'clone', 'size']
	fig = plt.figure(figsize=(10, 6))
	sns.scatterplot(data=dfu, x='clone', y='compartment', size=df.values.flatten(), alpha=0.8, sizes=(5, 4000))
	ax = fig.gca()
	ax.set_xlabel('Traced lineages')
	ax.set_ylabel('')
	ax.set_ylim(-0.6, len(df.columns) - 0.4)
	plt.xticks(rotation=90)
	plt.legend(bbox_to_anchor=(1.1, 0.55), loc=6, labelspacing=7, borderaxespad=0., scatterpoints=1, handletextpad=4, frameon=False)
	fig.savefig(args.pdf, bbox_inches='tight')
	logger.info('File %r written', args.pdf)


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
