"""
Some utility functions that work on sequences and lists of sequences.
"""
import random
import hashlib
import os
from collections import OrderedDict
import re
import numpy as np
import yaml
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


def sequence_hash(s, digits=4):
	"""
	For a string, return a 'fingerprint' that looks like 'S1234' (the character
	'S' is fixed). The idea is that this allows one to quickly see whether two
	sequences are not identical.
	"""
	h = int(hashlib.md5(s.encode()).hexdigest()[-4:], base=16)
	return 'S' + str(h % 10**digits).rjust(digits, '0')


class SerialPool:
	"""
	An alternative to multiprocessing.Pool that runs things in serial for
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
	"""Translate a nucleotide sequence to an amino acid sequence"""
	return ''.join(GENETIC_CODE.get(s[i:i+3], '*') for i in range(0, len(s), 3))


class UniqueNamer:
	"""
	Assign unique names by appending letters to already seen names.
	"""
	def __init__(self):
		self._names = set()

	def __call__(self, name):
		ext = 'A'
		new_name = name
		while new_name in self._names:
			if ext == '[':
				raise ValueError('Too many duplicate names')
			new_name = name + ext
			ext = chr(ord(ext) + 1)
		self._names.add(new_name)
		return new_name


class Merger:
	"""
	Merge similar items into one. To specify what "similar" means, implement
	the merged() method in a subclass.
	"""
	def __init__(self):
		self._items = []

	def add(self, item):
		# This method could possibly be made simpler if the graph structure
		# was made explicit.

		items = []
		m = None
		for existing_item in self._items:
			m = self.merged(existing_item, item)
			if m is None:
				items.append(existing_item)
			else:
				item = m
		items.append(item)
		self._items = items

	def __iter__(self):
		yield from self._items

	def merged(self, existing_item, item):
		"""
		If existing_item and item can be returned, this method must return
		a new item that represents a merged version of both. If they cannot
		be merged, it must return None.
		"""
		raise NotImplementedError("not implemented")


def relative_symlink(src, dst, force=False):
	"""
	Create a symbolic link in any directory.

	force -- if True, then overwrite an existing file/symlink
	"""
	if force:
		try:
			os.remove(dst)
		except FileNotFoundError:
			pass
	target = os.path.relpath(os.path.abspath(src), start=os.path.dirname(dst))
	os.symlink(target, dst)


class Config:
	PATH = 'igdiscover.yaml'
	def __init__(self, path=PATH):
		# Set some defaults.
		self.merge_program = 'pear'
		self.flash_maximum_overlap = 300
		self.limit = None  # or an integer
		self.cluster_program = 'vsearch'
		self.multialign_program = 'muscle-fast'
		self.maximum_expected_errors = None  # or an integer
		self.minimum_cluster_size = 200
		self.minimum_merged_read_length = 300
		self.mismatch_penalty = None
		self.barcode_length = 0
		self.iterations = 1
		self.ignore_j = False
		self.subsample = 500
		self.stranded = False
		self.forward_primers = None
		self.reverse_primers = None
		self.library_name = os.path.basename(os.getcwd())
		self.read_from(path)

	def read_from(self, path):
		with open(path) as f:
			content = f.read()
		self.__dict__.update(yaml.safe_load(content))
