"""
Run IgBLAST
"""
import sys
import os
import multiprocessing
import subprocess
from itertools import islice

from sqt import SequenceReader


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('igblast', help=__doc__)
	subparser.set_defaults(func=igblast_command)
	add = subparser.add_argument
	add('--threads', '-t', type=int, default=1,
		help='Number of threads (default: %(default)s)')
	add('--penalty', type=int, choices=(-1, -2, -3, -4), default=None,
		help='BLAST mismatch penalty (default: -1)')
	add('--species', default='rhesus_monkey',
		help='Which species (default: %(default)s)')
	add('--limit', type=int, metavar='N',
		help='Limit processing to first N records')
	add('database', help='path to database')
	add('fasta', help='File with original reads')
	return subparser


def igblast_command(args):
	"""
	Run IgBLAST in parallel
	"""
	if not 'IGDATA' in os.environ:
		raise ValueError("The IGDATA environment variable needs to be set")

	chunks = chunked_fasta(args.fasta, limit=args.limit)
	runner = Runner(args.database, args.species, args.penalty)
	with multiprocessing.Pool(args.threads) as pool:
		for result in pool.imap(runner, chunks, chunksize=1):
			sys.stdout.write(result)


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
	arguments += [
		#TODO '-auxiliary_data', '$IGDATA/optional_file/{species}_gl.aux',
		#'-auxiliary_data', '/dev/null',
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
