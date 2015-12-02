"""
Some utility functions that work on sequences and lists of sequences.
"""
import random
import hashlib
from collections import OrderedDict
import numpy as np
import re
from sqt.align import edit_distance, multialign, consensus
from sqt.dna import GENETIC_CODE


def downsampled(population, size):
	"""
	Return a random subsample of the population.

	Uses reservoir sampling. See https://en.wikipedia.org/wiki/Reservoir_sampling
	"""
	sample = population[:size]
	for index in range(size, len(population)):
		r = random.randint(0, index)
		if r < size:
			sample[r] = population[index]
	return sample


def distances(sequences, band=0.2):
	"""
	Compute all pairwise edit distances and return a square matrix.

	Entry [i,j] in the matrix is the edit distance between sequences[i]
	and sequences[j].
	"""
	m = np.zeros((len(sequences), len(sequences)), dtype=float)
	maxdiff = max((int(len(s) * band) for s in sequences), default=0)
	for i, s in enumerate(sequences):
		for j, t in enumerate(sequences):
			if i < j:
				m[j,i] = m[i,j] = min(maxdiff+1, edit_distance(s, t, maxdiff=maxdiff))
	return m


def iterative_consensus(sequences, program='muscle-medium', threshold=0.6,
		subsample_size=200, maximum_subsample_size=1600):
	"""
	Compute a consensus sequence of the given sequences, but do not use all
	sequences if there are many: First, try to compute the consensus from a
	small subsample. If there are 'N' bases, increase the subsample size and
	repeat until either there are no more 'N' bases, all available sequences
	have been used or maximum_subsample_size is reached.
	"""
	while True:
		sample = downsampled(sequences, subsample_size)
		aligned = multialign(OrderedDict(enumerate(sample)), program=program)
		cons = consensus(aligned, threshold=threshold).strip('N')
		if 'N' not in cons:
			# This consensus is good enough
			break
		if len(sequences) <= subsample_size:
			# We have already used all the sequences that are available
			break
		subsample_size *= 2
		if subsample_size > maximum_subsample_size:
			break
	return cons


def sequence_hash(s):
	"""
	Return a hash of a string that looks like 'S123' (S is fixed). The idea is
	that this allows one to quickly see whether two sequences are identical.
	"""
	h = int(hashlib.md5(s.encode()).hexdigest()[-6:], base=16)
	return 'S{:06}'.format(h % 1000000)


class SerialPool:
	"""
	An alternative to multiprocessing.Pool that runs things in parallel for
	easier debugging
	"""
	def __init__(self, *args, **kwargs):
		pass

	def __enter__(self):
		return self

	def __exit__(self, *args):
		pass

	def imap(self, func, iterable, chunksize):
		for i in iterable:
			yield func(i)


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
	"""
	Use this function as sorting key to sort in 'natural' order.

	>>> names = ['file10', 'file1.5', 'file1']
	>>> sorted(names, key=natural_sort_key)
	['file1', 'file1.5', 'file10']

	Source: http://stackoverflow.com/a/16090640/715090
	"""
	return [int(text) if text.isdigit() else text.lower()
		for text in re.split(_nsre, s)]


def nt_to_aa(s):
	"""Translate nucleotide sequence to amino acid sequence"""
	return ''.join(GENETIC_CODE.get(s[i:i+3], '*') for i in range(0, len(s), 3))
