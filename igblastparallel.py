#!/usr/bin/env python3
"""
Run IgBLAST in parallel.
"""
import sys
import subprocess
import multiprocessing
import logging
from sqt import HelpfulArgumentParser, SequenceReader

logger = logging.getLogger(__name__)

def run_igblast(fasta, database, organism='rhesus_monkey'):
	"""
	fasta -- path to a FASTA file or data in FASTA format as a string
	database -- directory that contains IgBLAST databases. Files in that
	directory must be databases created by the makeblastdb program and have
	names organism_gene, such as "rhesus_monkey_V".

	Only one of fasta_data and fasta_path may be given.

	Return IgBLAST’s output as a string.
	"""
	"""
	TODO
	Igblastwrapper has the databases in files named like this:
	$IGDATA/database/rhesus_monkey_IG_H_J
	Should we also include the _IG_H part?

	TODO
	When running igblastn, the IGDATA environment variable needs to point to
	the directory that contains the internal_data/ directory.
	Either require IGDATA to be set or set it here.
	"""
	if fasta.startswith('>'):
		fasta_input = '-'
	else:
		fasta_input = fasta
	arguments = ['igblastn']
	for gene in 'V', 'D', 'J':
		arguments += ['-germline_db_{gene}'.format(gene=gene),
			'{database}/{organism}_{gene}'.format(database=database, organism=organism, gene=gene)]
	arguments += [
		#TODO '-auxiliary_data', '$IGDATA/optional_file/{organism}_gl.aux',
		'-organism', organism,
		'-ig_seqtype', 'Ig',
		'-num_threads' , '1',
		'-domain_system', 'imgt',
		'-num_alignments_V', '1',
		'-num_alignments_D', '1',
		'-num_alignments_J', '1',
		'-outfmt', "7 qseqid qstart qseq sstart sseq pident",
		'-query', fasta_input,
		'-out', '-',  # write to stdout
	]
	result = subprocess.check_output(arguments, input=fasta if fasta_input == '-' else None, universal_newlines=True)
	return result


def chunked_fasta(path, chunksize=500):
	"""
	Split a FASTA file into chunks that contain at most chunksize records.
	"""
	def buf_to_fasta(b):
		return ''.join(">{}\n{}\n".format(r.name, r.sequence) for r in b)
	sr = SequenceReader(path)
	buf = []
	for record in sr:
		buf.append(record)
		if len(buf) == chunksize:
			yield buf_to_fasta(buf)
			buf = []
	if len(buf) > 0:
		yield buf_to_fasta(buf)
	sr.close()


class Runner:
	"""nested functions cannot be pickled"""
	def __init__(self, dbpath):
		self.dbpath = dbpath

	def __call__(self, s):
		return run_igblast(s, self.dbpath)


def main():
	logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('--threads', type=int, default=1, help='Number of threads')
	parser.add_argument('database', help='path to database')
	parser.add_argument('fasta', help='FASTA file with sequences')
	args = parser.parse_args()

	def my_run(s):
		return run_igblast(s, 'pipeline/database')
	#sys.stdout.write(run_igblast(sys.argv[1], 'pipeline/database'))
	chunks = chunked_fasta(args.fasta)
	runner = Runner(args.database)
	with multiprocessing.Pool(args.threads) as pool:
		for result in pool.map(runner, chunks):
			sys.stdout.write(result)


if __name__ == '__main__':
	main()
