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
from sqt.dna import GENETIC_CODE, nt_to_aa as _nt_to_aa


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
		yield from sorted(self._items, key=lambda x: x.name)

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


class ConfigurationError(Exception):
	pass


class Config:
	DEFAULT_PATH = 'igdiscover.yaml'

	def __init__(self, file):
		# Set some defaults.
		self.debug = False
		self.species = None
		self.merge_program = 'pear'
		self.flash_maximum_overlap = 300
		self.limit = None  # or an integer
		self.cluster_program = 'vsearch'
		self.multialign_program = 'muscle-fast'
		self.maximum_expected_errors = None  # or an integer
		self.minimum_merged_read_length = 300
		self.mismatch_penalty = None
		self.barcode_length = 0
		self.iterations = 3
		self.ignore_j = False
		self.subsample = 1000
		self.stranded = False
		self.forward_primers = None
		self.reverse_primers = None
		self.rename = True
		self.race_g = False
		self.seed = 1
		self.exact_copies = None
		self.preprocessing_filter = dict(v_coverage=90, j_coverage=60, v_evalue=1E-3)
		self.pre_germline_filter = dict(unique_cdr3s=2, unique_js=2, check_motifs=False,
			whitelist=True, cluster_size=0, differences=1, allow_stop=True, cross_mapping_ratio=0.02)
		self.germline_filter = dict(unique_cdr3s=5, unique_js=3, check_motifs=False,
			whitelist=True, cluster_size=100, differences=1, allow_stop=False, cross_mapping_ratio=0.02)
		self.cdr3_location = [-80, -60]
		self.library_name = os.path.basename(os.getcwd())

		self.read_from(file)

	def read_from(self, file):
		content = file.read()
		new_config = self.make_compatible(yaml.safe_load(content))
		for key in ('preprocessing_filter', 'pre_germline_filter', 'germline_filter'):
			if key in new_config:
				self.__dict__[key].update(new_config[key])
				del new_config[key]
		self.__dict__.update(new_config)

	def make_compatible(self, config):
		"""
		Convert old-style configuration to new style. Raise ConfigurationError if configuration is invalid.
		Return updated config dict.
		"""
		if 'barcode_length' in config and ('barcode_length_5prime' in config or 'barcode_length_3prime' in config):
			raise ConfigurationError('Old-style configuration of barcode length via "barcode_length" option cannot be '
					'used at the same time as new-style configuration via "barcode_length_5prime" or "barcode_length_3prime"')
		barcode_length_5prime = config.get('barcode_length_5prime', 0)
		barcode_length_3prime = config.get('barcode_length_3prime', 0)
		if barcode_length_5prime > 0 and barcode_length_3prime > 0:
			raise ConfigurationError('barcode_length_5prime and barcode_length_3prime can currently not both be greater than zero.')
		if barcode_length_5prime > 0:
			config['barcode_length'] = barcode_length_5prime
		elif barcode_length_3prime > 0:
			config['barcode_length'] = -barcode_length_3prime
		config.pop('barcode_length_5prime', None)
		config.pop('barcode_length_3prime', None)

		if 'seed' in config and config['seed'] is False:
			config['seed'] = None
		return config

	@classmethod
	def from_default_path(cls):
		with open(cls.DEFAULT_PATH) as f:
			return Config(file=f)


def nt_to_aa(s):
	"""Translate a nucleotide sequence to an amino acid sequence"""
	try:
		# try fast version first
		return _nt_to_aa(s)
	except ValueError:
		# failure because there was an unknown nucleotide
		return ''.join(GENETIC_CODE.get(s[i:i+3], '*') for i in range(0, len(s), 3))


def has_stop(sequence):
	"""
	Return a boolean indicating whether the sequence has an internal stop codon.
	An incomplete codon at the end is allowed.

	>>> has_stop('GGG')
	False
	>>> has_stop('TAA')
	True
	>>> has_stop('GGGAC')
	False
	"""
	s = sequence[:len(sequence) // 3 * 3]
	return '*' in nt_to_aa(s)


def plural_s(n):
	return 's' if n != 1 else ''
