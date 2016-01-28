"""
Rename sequences in a target FASTA file using a template FASTA file

The sequences in the target file will get the name that they have in the
template file. Sequences are considered to be equivalent if one is a prefix of
the other.
"""
import logging
from collections import Counter
from sqt import FastaReader
from .utils import natural_sort_key

logger = logging.getLogger(__name__)

def add_arguments(parser):
	arg = parser.add_argument
	arg('--no-sort', dest='sort', action='store_false', default=True,
		help='Do not sort sequences by name')
	arg('--not-found', metavar='TEXT', default=' (not found)',
		help='Append this text to the record name when the sequence was not found '
		'in the template. Default: %(default)r')
	arg('template', help='FASTA file with correctly named sequences')
	arg('target', help='FASTA file with to-be renamed sequences')


class PrefixDict:
	"""
	A dict that maps strings to values, but where a prefix of a key is enough
	to retrieve the value.
	"""
	def __init__(self, items):
		self._items = []
		for k, v in items:
			self.add(k, v)

	def add(self, k, v):
		try:
			self[k]
			raise ValueError('Key {!r} already exists'.format(k))
		except KeyError:
			self._items.append((k, v))

	def __getitem__(self, key):
		found = None
		for seq, value in self._items:
			if seq.startswith(key) or key.startswith(seq):
				if found is not None:
					# TODO don't use keyerror here
					raise KeyError('Key {!r} is ambiguous'.format(key))
				found = value
		if found is None:
			raise KeyError(key)
		return found

	def get(self, key, default=None):
		try:
			v = self[key]
		except KeyError:
			return default
		else:
			return v

	def __len__(self):
		return len(self._items)

def main(args):
	with FastaReader(args.template) as fr:
		template = PrefixDict([])
		for record in fr:
			try:
				template.add(record.sequence.upper(), record.name)
			except ValueError:
				logger.error('Sequences in entry %r and %r are duplicate',
					record.name, template[record.sequence.upper()])
	logger.info('Read %d entries from template', len(template))

	with FastaReader(args.target) as fr:
		sequences = list(fr)

	# Rename
	for record in sequences:
		record.name = template.get(record.sequence.upper(), record.name + args.not_found)

	if args.sort:
		sequences = sorted(sequences, key=lambda s: natural_sort_key(s.name))
	for record in sequences:
		print('>{}\n{}'.format(record.name, record.sequence))
	logger.info('Wrote %s FASTA records', len(sequences))
