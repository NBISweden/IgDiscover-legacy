"""
Draw a dendrogram of sequences in a FASTA file.
"""
import logging
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sqt import FastaReader
from igypipe.utils import distances


logger = logging.getLogger(__name__)

def add_arguments(parser):
	parser.add_argument('fasta', help='Path to input FASTA file')
	parser.add_argument('plot', help='Path to output PDF or PNG')


def main(args):
	sequences = list(FastaReader(args.fasta))
	m = distances([s.sequence for s in sequences])
	y = distance.squareform(m)
	l = hierarchy.average(y)  # UPGMA
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
	result = hierarchy.dendrogram(l, labels=labels, leaf_font_size=font_size, orientation='right', color_threshold=0.95*max(l[:,2]))
	ax.grid(False)
	fig.tight_layout()
	fig.savefig(args.plot)
