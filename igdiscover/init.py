"""
Create and initialize a new analysis directory.
"""
import glob
import logging
import re
import os
import os.path
import sys
import subprocess
import pkg_resources
from sqt import SequenceReader
from .config import Config

try:
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
except ImportError:
	tk = None

from xopen import xopen

logger = logging.getLogger(__name__)

do_not_show_cpustats = 1


def add_arguments(parser):
	parser.add_argument('--database',  metavar='PATH', default=None,
		help='Directory with V.fasta, D.fasta and J.fasta files. If not given, a dialog is shown.')
	parser.add_argument('--species', '--sp', metavar='SPECIES', default='any', choices=('any', 'human'),
		help='Use configuration settings specific to %(choices)s. Default: %(default)s.')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('--single-reads', default=None, metavar='READS',
		help='File with single-end reads (.fasta.gz or .fastq.gz)')
	group.add_argument('--reads1', default=None,
		help='First paired-end read file. The second is found automatically. '
			'Must be a .fastq.gz file. If not given, a dialog is shown.')
	parser.add_argument('directory', help='New analysis directory to create')


def launch(path):
	if hasattr(os, 'startfile'):
		os.startfile(path)
	elif sys.platform == 'linux':
		subprocess.call(['xdg-open', path])
	elif sys.platform == 'darwin':
		subprocess.call(['open', path])


class TkinterGui:
	"""Show a GUI for selecting reads and the database directory"""
	def __init__(self):
		if not tk:  # import failed
			raise ImportError()
		root = tk.Tk()
		root.withdraw()

	def yesno(self, title, question):
		return tk.messagebox.askyesno(title, question)

	def database_path(self, initialdir):
		path = tk.filedialog.askdirectory(
			title="Choose V/D/J database directory", mustexist=True,
			initialdir=initialdir)
		return path

	def reads1_path(self):
		path = tk.filedialog.askopenfilename(
			title="Choose first reads file",
			filetypes=[
				("Reads", "*.fastq.gz"),
				("Any file", "*")])
		return path

	def single_reads_path(self):
		path = tk.filedialog.askopenfilename(
			title="Choose single-end reads file",
			filetypes=[
				("Reads", "*.fasta *.fastq *.fastq.gz *.fasta.gz"),
				("Any file", "*")])
		return path

	def error(self, title, message):
		tk.messagebox.showerror(title, message)


"""
# Works, but let’s not introduce the PySide dependency for now.

def qt_path():
	import PySide
	from PySide.QtGui import QApplication
	from PySide.QtGui import QMessageBox, QFileDialog

	# Create the application object
	app = QApplication([])

	path = QFileDialog.getOpenFileName(None,
		"Open first reads file", '.', "FASTA/FASTQ reads (*.fastq *.fasta *.fastq.gz *.fasta.gz);; Any file (*)")
	# QMessageBox.information(None, 'Chosen file', path[0])
	return path[0]
"""


def is_1_2(s, t):
	"""
	Determine whether s and t are identical except for a single character of
	which one of them is '1' and the other is '2'.
	"""
	differences = 0
	one_two = {'1', '2'}
	for c1, c2 in zip(s, t):
		if c1 != c2:
			differences += 1
			if differences == 2:
				return False
			if {c1, c2} != one_two:
				return False
	return differences == 1


def guess_paired_path(path):
	"""
	Given the path to a file that contains the sequences for the first read in a
	pair, return the file that contains the sequences for the second read in a
	pair. Both files must have identical names, except that the first must have
	a '1' in its name, and the second must have a '2' at the same position.

	Return None if no second file was found or if there are too many candidates.

	>>> guess_paired_path('file.1.fastq.gz')  # doctest: +SKIP
	'file.2.fastq.gz'  # if that file exists
	"""
	base, name = os.path.split(path)
	# All lone 1 digits replaced with '?'
	name_with_globs = re.sub(r'(?<![0-9])1(?![0-9])', '?', name)
	glob_pattern = os.path.join(base, name_with_globs)
	paths = [p for p in glob.glob(glob_pattern) if is_1_2(p, path)]
	if len(paths) == 1:
		return paths[0]
	return None


class UnknownFileFormatError(Exception):
	pass


def file_type(path):
	"""
	Return 'fasta' or 'fastq' depending on file format. The file may optionally
	be gzip-compressed.
	"""
	if path.endswith('.gz'):
		file = xopen(path)
	else:
		file = open(path)
	with file as f:
		first_char = f.read(1)
		if first_char == '@':
			return 'fastq'
		elif first_char == '>':
			return 'fasta'
		else:
			raise UnknownFileFormatError('Cannot recognize format. File starts with neither ">" nor "@"')


def try_open(path):
	try:
		with open(path) as f:
			pass
	except OSError as e:
		logger.error('Could not open %r: %s', path, e)
		sys.exit(1)


def read_and_repair_fasta(path):
	"""
	Read a FASTA file and make sure it is suitable for use with makeblastdb.
	It repairs the following issues:

	- If a record is empty, it is skipped
	- If a record name occurs more than once, the second record name gets a suffix
	- If a sequence occurs more than once, occurrences after the first are skipped
	"""
	with SequenceReader(path) as sr:
		records = list(sr)

	names = set()
	sequences = dict()
	for r in records:
		r.sequence = r.sequence.upper()
		if len(r.sequence) == 0:
			logger.info("Record %r is empty, skipping it.", r.name)
			continue
		name = r.name
		i = 0
		while name in names:
			i += 1
			name = r.name + '_{}'.format(i)
		if name != r.name:
			logger.info('Record name %r replaced with %r because it occurs more than once', r.name, name)
		if r.sequence in sequences:
			logger.info('Skipping %r because it contains the same sequence as %r',
				r.name, sequences[r.sequence])
			continue
		sequences[r.sequence] = name
		names.add(name)
		r.name = name
		yield r


def main(args):
	if ' ' in args.directory:
		logger.error('The name of the analysis directory must not contain spaces')
		sys.exit(1)

	if os.path.exists(args.directory):
		logger.error('The target directory {!r} already exists.'.format(args.directory))
		sys.exit(1)

	# If reads files or database were not given, initialize the GUI
	if (args.reads1 is None and args.single_reads is None) or args.database is None:
		try:
			gui = TkinterGui()
		except ImportError:  # TODO tk.TclError cannot be caught when import of tk fails
			logger.error('GUI cannot be started. Please provide reads1 file '
				'and database directory on command line.')
			sys.exit(1)
	else:
		gui = None

	# Find out whether data is paired or single
	assert not (args.reads1 and args.single_reads)
	if args.reads1 is None and args.single_reads is None:
		paired = gui.yesno('Paired end or single-end reads',
			'Are your reads paired and need to be merged?\n\n'
			'If you answer "Yes", next select the FASTQ files '
			'with the first of your paired-end reads.\n'
			'If you answer "No", next select the FASTA or FASTQ '
			'file with single-end reads.')
		if paired is None:
			logger.error('Cancelled')
			sys.exit(2)
	else:
		paired = bool(args.reads1)

	# Assign reads1 and (if paired) also reads2
	if paired:
		if args.reads1 is not None:
			reads1 = args.reads1
			try_open(reads1)
		else:
			reads1 = gui.reads1_path()
			if not reads1:
				logger.error('Cancelled')
				sys.exit(2)
		reads2 = guess_paired_path(reads1)
		if reads2 is None:
			logger.error('Could not determine second file of paired-end reads')
			sys.exit(1)
	else:
		if args.single_reads is not None:
			reads1 = args.single_reads
			try_open(reads1)
		else:
			reads1 = gui.single_reads_path()
			if not reads1:
				logger.error('Cancelled')
				sys.exit(2)

	if args.database is not None:
		dbpath = args.database
	else:
		# TODO as soon as we distribute our own database files, we can use this:
		# database_path = pkg_resources.resource_filename('igdiscover', 'databases')
		databases_path = None
		dbpath = gui.database_path(databases_path)
		if not dbpath:
			logger.error('Cancelled')
			sys.exit(2)

	database = dict()
	for g in ['V', 'D', 'J']:
		path = os.path.join(dbpath, g + '.fasta')
		if not os.path.exists(path):
			logger.error(
				'The database directory %r must contain the three files '
				'V.fasta, D.fasta and J.fasta', dbpath)
			logger.error(
				'A dummy D.fasta is necessary even if analyzing light chains (see manual)')
			sys.exit(2)
		database[g] = list(read_and_repair_fasta(path))

	# Create the directory
	try:
		os.mkdir(args.directory)
	except OSError as e:
		logger.error(e)
		sys.exit(1)

	def create_symlink(readspath, dirname, target):
		gz = '.gz' if readspath.endswith('.gz') else ''
		if not os.path.isabs(readspath):
			src = os.path.relpath(readspath, dirname)
		else:
			src = readspath
		os.symlink(src, os.path.join(dirname, target + gz))

	if paired:
		create_symlink(reads1, args.directory, 'reads.1.fastq')
		create_symlink(reads2, args.directory, 'reads.2.fastq')
	else:
		try:
			target = 'reads.' + file_type(reads1)
		except UnknownFileFormatError:
			logger.error('Cannot determine whether reads file is FASTA or FASTQ')
			sys.exit(1)
		create_symlink(reads1, args.directory, target)

	configuration = pkg_resources.resource_string('igdiscover', Config.DEFAULT_PATH).decode()
	# Write the configuration file
	local_config_path = os.path.join(args.directory, Config.DEFAULT_PATH)
	with open(local_config_path, 'w') as f:
		f.write(configuration)

	# Update defaults that are human specific
	if args.species.lower() == 'human':
		Config.update_option(local_config_path, [('preprocessing_filter.v_coverage', '97')])

	# Create database files
	database_dir = os.path.join(args.directory, 'database')
	os.mkdir(database_dir)
	for gene in ['V', 'D', 'J']:
		with open(os.path.join(database_dir, gene + '.fasta'), 'w') as db_file:
			for record in database[gene]:
				print('>{}\n{}'.format(record.name, record.sequence), file=db_file)

	if gui is not None:
		# Only suggest to edit the config file if at least one GUI dialog has been shown
		if gui.yesno('Directory initialized',
				'Do you want to edit the configuration file now?'):
			launch(os.path.join(args.directory, Config.DEFAULT_PATH))
	logger.info('Directory %s initialized.', args.directory)
	logger.info('Edit %s/%s, then run "cd %s && igdiscover run" to start the analysis', args.directory, Config.DEFAULT_PATH, args.directory)
