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

from sqt import SequenceReader
from sqt.utils import available_cpu_count


def add_arguments(parser):
	add = parser.add_argument
	add('--threads', '-t', '-j', type=int, default=available_cpu_count(),
		help='Number of threads. Default: no. of available CPUs (%(default)s)')
	add('--penalty', type=int, choices=(-1, -2, -3, -4), default=None,
		help='BLAST mismatch penalty (default: -1)')
	add('--species', default='rhesus_monkey',
		help='Which species (default: %(default)s)')
	add('--limit', type=int, metavar='N',
		help='Limit processing to first N records')
	add('database', help='Database directory. Must contain BLAST databases for '
		'SPECIES_V, SPECIES_D and SPECIES_J.')
	add('fasta', help='File with original reads')


def run_igblast(fasta, database, species, penalty=None):
	"""
	fasta -- path to a FASTA file or data in FASTA format as a string
	database -- directory that contains IgBLAST databases. Files in that
	directory must be databases created by the makeblastdb program and have
	names organism_gene, such as "rhesus_monkey_V".

	Return IgBLAST’s output as a string.
	"""

	"""
	TODO
	Igblastwrapper has the databases in files named like this:
	$IGDATA/database/rhesus_monkey_IG_H_J
	Should we also include the _IG_H part?
	"""
	if fasta.startswith('>'):
		fasta_input = '-'
	else:
		fasta_input = fasta
	arguments = ['igblastn']
	for gene in 'V', 'D', 'J':
		arguments += ['-germline_db_{gene}'.format(gene=gene),
			os.path.join(database, '{species}_{gene}'.format(species=species, gene=gene))]
	if penalty is not None:
		arguments += ['-penalty', str(penalty)]
	# The empty aux suppresses a warning from IgBLAST. /dev/null does not work.
	empty_aux_path = pkg_resources.resource_filename('igdiscover', 'empty.aux')
	arguments += [
		'-auxiliary_data', empty_aux_path,
		'-organism', species,
		'-ig_seqtype', 'Ig',
		'-num_threads' , '1',
		'-domain_system', 'imgt',
		'-num_alignments_V', '1',
		'-num_alignments_D', '1',
		'-num_alignments_J', '1',
		'-outfmt', "7 sseqid qstart qseq sstart sseq pident slen evalue",
		'-query', fasta_input,
		'-out', '-',  # write to stdout
	]
	result = subprocess.check_output(arguments, input=fasta if fasta_input == '-' else None, universal_newlines=True)
	return result


def chunked_fasta(path, chunksize=500, limit=None):
	"""
	Split a FASTA file into chunks that contain at most chunksize records.
	Take only the first 'limit' records, or all if limit is None.
	"""
	def buf_to_fasta(b):
		return ''.join(">{}\n{}\n".format(r.name, r.sequence) for r in b)
	sr = SequenceReader(path)
	buf = []
	for record in islice(sr, 0, limit):
		buf.append(record)
		if len(buf) == chunksize:
			yield buf_to_fasta(buf)
			buf = []
	if len(buf) > 0:
		yield buf_to_fasta(buf)
	sr.close()


class Runner:
	"""nested functions cannot be pickled"""
	def __init__(self, dbpath, species, penalty):
		self.dbpath = dbpath
		self.species = species
		self.penalty = penalty

	def __call__(self, s):
		return run_igblast(s, self.dbpath, self.species, self.penalty)


def main(args):
	chunks = chunked_fasta(args.fasta, limit=args.limit)
	runner = Runner(args.database, args.species, args.penalty)
	with multiprocessing.Pool(args.threads) as pool:
		for result in pool.imap(runner, chunks, chunksize=1):
			sys.stdout.write(result)
