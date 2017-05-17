"""
Run IgBLAST and output a result table

This is a wrapper for the "igblastn" tool which has a simpler command-line
syntax and can also run IgBLAST in parallel.

The results are parsed, postprocessed and printed as a tab-separated table
to standard output.

Postprocessing includes:

- The CDR3 is detected by using a regular expression
- The leader is detected within the sequence before the found V gene (by
  searching for the start codon).
- If the V sequence hit starts not at base 1 in the reference, it is extended
  to the left.
"""
import sys
import os
import multiprocessing
import subprocess
from contextlib import ExitStack
from io import StringIO
from itertools import islice

import errno
import pkg_resources
import logging
import tempfile
import json

from xopen import xopen

from sqt import SequenceReader
from sqt.utils import available_cpu_count

from .utils import get_cpu_time, SerialPool
from .parse import TableWriter, ExtendedIgBlastParser


logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--threads', '-t', '-j', type=int, default=available_cpu_count(),
		help='Number of threads. Default: no. of available CPUs (%(default)s)')
	arg('--penalty', type=int, choices=(-1, -2, -3, -4), default=None,
		help='BLAST mismatch penalty (default: -1)')
	arg('--species', default='mouse',
		help='Which species. Note that this setting does not seem to have '
			'any effect since we provide our own database to IgBLAST. '
			'Default: %(default)s')
	arg('--raw', metavar='FILE', help='Write raw IgBLAST output to FILE '
			'(add .gz to compress)')
	arg('--limit', type=int, metavar='N',
		help='Limit processing to first N records')
	arg('--rename', default=None, metavar='PREFIX',
		help='Rename reads to PREFIXseqN (where N is a number starting at 1)')
	arg('--stats', metavar='FILE',
		help='Write statistics in JSON format to FILE')

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


def chunked(iterable, chunksize: int):
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

	It runs IgBLAST and parses the output for a list of sequences.
	"""
	def __init__(self, dbpath, species, penalty, v_database):
		self.dbpath = dbpath
		self.species = species
		self.penalty = penalty
		self.v_database = v_database

	def __call__(self, sequences):
		"""
		Return tuples (igblast_result, records) where igblast_result is the raw IgBLAST output
		and records is a list of ExtendedIgBlastRecord objects (the parsed output).
		"""
		igblast_result = run_igblast(sequences, self.dbpath, self.species, self.penalty)
		parser = ExtendedIgBlastParser(sequences,
			StringIO(igblast_result), self.v_database)
		records = list(parser)
		assert len(records) == len(sequences)
		return igblast_result, records


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
	with SequenceReader(os.path.join(args.database, 'V.fasta')) as sr:
		v_database = {record.name: record.sequence.upper() for record in sr}

	detected_cdr3s = 0
	writer = TableWriter(sys.stdout)
	with ExitStack() as stack:
		if args.raw:
			raw_output = stack.enter_context(xopen(args.raw, 'w'))
		else:
			raw_output = None
		tmpdir = stack.enter_context(tempfile.TemporaryDirectory())

		# Create the three BLAST databases in a temporary directory
		for gene in list('VDJ'):
			makeblastdb(os.path.join(args.database, gene + '.fasta'), os.path.join(tmpdir, gene))

		fasta = stack.enter_context(SequenceReader(args.fasta))
		chunks = chunked(islice(fasta, 0, args.limit), chunksize=1000)
		runner = Runner(tmpdir, args.species, args.penalty, v_database)
		pool = stack.enter_context(multiprocessing.Pool(args.threads))
		n = 0  # number of records processed so far
		for igblast_output, igblast_records in pool.imap(runner, chunks, chunksize=1):
			if raw_output:
				raw_output.write(igblast_output)
			for record in igblast_records:
				n += 1
				if args.rename is not None:
					record.query_name = "{}seq{}".format(args.rename, n)
				d = record.asdict()
				if d['CDR3_aa']:
					detected_cdr3s += 1
				try:
					writer.write(d)
				except IOError as e:
					if e.errno == errno.EPIPE:
						sys.exit(1)
					raise

	cpu_time = get_cpu_time()
	if cpu_time is not None:
		logger.info('Processed {} sequences at {:.1f} ms/sequence'.format(n, cpu_time / n * 1E3))

	logger.info('%d IgBLAST assignments parsed and written', n)
	logger.info('CDR3s detected in %.1f%% of all sequences', detected_cdr3s / n * 100)
	if args.stats:
		stats = {'total': n, 'detected_cdr3s': detected_cdr3s}
		with open(args.stats, 'w') as f:
			json.dump(stats, f)
			print(file=f)
