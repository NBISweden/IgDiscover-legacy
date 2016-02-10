"""
Draw a dendrogram of sequences in a FASTA file.
"""
import logging
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sqt import FastaReader
from igdiscover.utils import distances


logger = logging.getLogger(__name__)

def add_arguments(parser):
	parser.add_argument('--mark', '--db',
		help='Path to a FASTA file with a set of "known" sequences. Sequences '
		'in the main file that do *not* occur here will be marked with (new). '
		'If not given, no sequences will be marked (use this to compare two '
		'databases.')
	parser.add_argument('fasta', help='Path to input FASTA file')
	parser.add_argument('plot', help='Path to output PDF or PNG')


class PrefixComparer:
	def __init__(self, sequences):
		self._sequences = [ s.upper() for s in sequences ]

	def __contains__(self, other):
		for seq in self._sequences:
			if seq.startswith(other) or other.startswith(seq):
				return True
		return False


def main(args):
	with FastaReader(args.fasta) as fr:
		sequences = list(fr)
	logger.info('Plotting dendrogram of %s sequences', len(sequences))
	if args.mark:
		with FastaReader(args.mark) as fr:
			mark = PrefixComparer(record.sequence for record in fr)
		labels = []
		n_new = 0
		for record in sequences:
			if record.sequence not in mark:
				extra = ' (new)'
				n_new += 1
			else:
				extra = ''
			labels.append(record.name + extra)
		logger.info('%s sequence(s) marked as "new"', n_new)
	else:
		labels = [s.name for s in sequences]
	sns.set_style("white")
	font_size = 297 / 25.4 * 72 / (len(labels) + 5)
	font_size = min(16, max(6, font_size))
	height = font_size * (len(labels) + 5) / 72
	fig = plt.figure(figsize=(210 / 25.4, height))
	matplotlib.rcParams.update({'font.size': 4})
	ax = fig.gca()
	sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
	sns.set_style('whitegrid')
	if len(sequences) >= 2:
		m = distances([s.sequence for s in sequences])
		y = distance.squareform(m)
		mindist = int(y.min())
		logger.info('Smallest distance is %s. Found between:', mindist)
		for i,j in np.argwhere(m == y.min()):
			if i < j:
				logger.info('%s and %s', labels[i], labels[j])
		l = hierarchy.average(y)  # UPGMA
		hierarchy.dendrogram(l, labels=labels, leaf_font_size=font_size, orientation='right', color_threshold=0.95*max(l[:,2]))
	else:
		ax.text(0.5, 0.5, 'no sequences', fontsize='xx-large')
	ax.grid(False)
	fig.set_tight_layout(True)
	fig.savefig(args.plot)
