"""
Run IgBLAST

This is a wrapper for the igblastn tool that has a somewhat simpler command-line
syntax and can also run IgBLAST in parallel. The result is printed to standard
output.
"""
import sys
import os
import multiprocessing
import subprocess
from itertools import islice
import pkg_resources
import logging
import tempfile

from sqt import SequenceReader
from sqt.utils import available_cpu_count

from .utils import get_cpu_time


logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--threads', '-t', '-j', type=int, default=available_cpu_count(),
		help='Number of threads. Default: no. of available CPUs (%(default)s)')
	arg('--penalty', type=int, choices=(-1, -2, -3, -4), default=None,
		help='BLAST mismatch penalty (default: -1)')
	arg('--species', default='mouse',
		help='Which species (default: %(default)s)')
	arg('--limit', type=int, metavar='N',
		help='Limit processing to first N records')
	arg('database', help='Database directory with V.fasta, D.fasta, J.fasta.')
	arg('fasta', help='File with original reads')


def run_igblast(sequences, database, species, penalty=None):
	"""
	sequences -- list of Sequence objects
	database -- directory that contains IgBLAST databases. Files in that
	directory must be databases created by the makeblastdb program and have
	names V, D, and J.

	Return IgBLAST’s raw output as a string.
	"""
	arguments = ['igblastn']
	for gene in 'V', 'D', 'J':
		arguments += ['-germline_db_{gene}'.format(gene=gene),
			os.path.join(database, '{gene}'.format(gene=gene))]
	if penalty is not None:
		arguments += ['-penalty', str(penalty)]
	# An empty .aux suppresses a warning from IgBLAST. /dev/null does not work.
	empty_aux_path = pkg_resources.resource_filename('igdiscover', 'empty.aux')
	arguments += [
		'-auxiliary_data', empty_aux_path,
		'-organism', species,
		'-ig_seqtype', 'Ig',
		'-num_threads', '1',
		'-domain_system', 'imgt',
		'-num_alignments_V', '1',
		'-num_alignments_D', '1',
		'-num_alignments_J', '1',
		'-outfmt', '7 sseqid qstart qseq sstart sseq pident slen evalue',
		'-query', '-',
		'-out', '-',  # write to stdout
	]
	fasta_str = ''.join(">{}\n{}\n".format(r.name, r.sequence) for r in sequences)
	return subprocess.check_output(arguments, input=fasta_str, universal_newlines=True)


def chunked(iterable, chunksize):
	"""
	Group the iterable into lists of length chunksize
	>>> list(chunked('ABCDEFG', 3))
	[['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
	"""
	chunk = []
	for it in iterable:
		if len(chunk) == chunksize:
			yield chunk
			chunk = []
		chunk.append(it)
	if chunk:
		yield chunk


class Runner:
	"""
	This is the target of a multiprocessing pool. The target needs to
	be pickleable, and because nested functions cannot be pickled,
	we need this separate class.
	"""
	def __init__(self, dbpath, species, penalty):
		self.dbpath = dbpath
		self.species = species
		self.penalty = penalty

	def __call__(self, chunk):
		return len(chunk), run_igblast(chunk, self.dbpath, self.species, self.penalty)


def makeblastdb(fasta, database_name):
	with SequenceReader(fasta) as fr:
		sequences = list(fr)
	if not sequences:
		raise ValueError("FASTA file {} is empty".format(fasta))
	process_output = subprocess.check_output(
		['makeblastdb', '-parse_seqids', '-dbtype', 'nucl', '-in', fasta, '-out', database_name],
		stderr=subprocess.STDOUT
	)
	if b'Error: ' in process_output:
		raise subprocess.SubprocessError()


def main(args):
	with tempfile.TemporaryDirectory() as tmpdir:
		# Create the three BLAST databases in a temporary directory
		for gene in list('VDJ'):
			makeblastdb(os.path.join(args.database, gene + '.fasta'), os.path.join(tmpdir, gene))

		with SequenceReader(args.fasta) as fasta:
			chunks = chunked(islice(fasta, 0, args.limit), chunksize=1000)
			runner = Runner(tmpdir, args.species, args.penalty)
			with multiprocessing.Pool(args.threads) as pool:
				total = 0
				for n, igblast_output in pool.imap(runner, chunks, chunksize=1):
					sys.stdout.write(igblast_output)
					total += n

	cpu_time = get_cpu_time()
	if cpu_time is not None:
		logger.info('Processed {} sequences at {:.1f} ms/sequence'.format(total, cpu_time / total * 1E3))
