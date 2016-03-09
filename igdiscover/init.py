"""
Create and initialize a new analysis directory.
"""
import glob
import logging
import os
import os.path
import sys
import shutil
import subprocess
import pkg_resources
try:
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
except ImportError:
	tk = None

from cutadapt.xopen import xopen

logger = logging.getLogger(__name__)


PIPELINE_CONF = 'igdiscover.yaml'


def add_arguments(parser):
	parser.add_argument('--database', '--db', metavar='PATH', default=None,
		help='Directory with IgBLAST database files. If not given, a dialog is shown.')
	group = parser.add_mutually_exclusive_group()
	group.add_argument('--single-reads', default=None, metavar='READS',
		help='File with single-end reads (.fasta.gz or .fastq.gz)')
	group.add_argument('--reads1', default=None,
		help='First paired-end read file. The second is found automatically. '
			'Must be a .fastq.gz file. If not given, a dialog is shown.')
	parser.add_argument('--library-name', metavar='NAME', default=None,
		help='Name of the library. Sets library_name in the configuration file.')
	parser.add_argument('directory', help='New analysis directory to create')


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
			if set([c1, c2]) != one_two:
				return False
	return differences == 1


def guess_paired_path(path):
	"""
	Given the path to a file that contains the sequences for the first read in a
	pair, return the file that contains the sequences for the second read in a
	pair. Both files must have identical names, except that the first must have
	a '1' in its name, and the second must have a '2' at the same position.

	Return None if no second file was found or if there are too many candidates.

	>>> guess_paired_path('file.1.fastq.gz')
	'file.2.fastq.gz'  # if that file exists
	"""
	base, name = os.path.split(path)
	glob_pattern = os.path.join(base, name.replace('1', '?'))
	paths = [ p for p in glob.glob(glob_pattern) if is_1_2(p, path) and '_R1_' not in p ]
	if len(paths) != 1:
		return None
	return paths[0]


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


def main(args):
	if ' ' in args.directory:
		sys.exit('The name of the analysis directory must not contain spaces')

	# If reads files or database were not given, initialize the GUI
	if (args.reads1 is None and args.single_reads is None) or args.database is None:
		try:
			gui = TkinterGui()
		except (ImportError, tk.TclError):
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
			'with the <em>first</em> of your paired-end reads.\n'
			'If you answer "No", next select the FASTA or FASTQ '
			'file with your sequences.')
		if paired is None:
			logger.error('Cancelled')
			sys.exit(2)
	else:
		paired = bool(args.reads1)

	# Assign reads1 and (if paired) also reads2
	if paired:
		if args.reads1 is not None:
			reads1 = args.reads1
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

	# Create the directory
	try:
		os.mkdir(args.directory)
	except OSError as e:
		logger.error(e)
		sys.exit(1)

	def create_symlink(readspath, dirname, target):
		gz = '.gz' if readspath.endswith('.gz') else ''
		rel = os.path.relpath(readspath, dirname)
		os.symlink(rel, os.path.join(dirname, target + gz))

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

	if args.library_name:
		library_name = args.library_name
	else:
		library_name = os.path.basename(os.path.normpath(args.directory))

	# Write the configuration file
	configuration = pkg_resources.resource_string('igdiscover', PIPELINE_CONF).decode()
	with open(os.path.join(args.directory, PIPELINE_CONF), 'w') as f:
		for line in configuration.splitlines(keepends=True):
			if line.startswith('library_name:'):
				line = 'library_name: ' + library_name + '\n'
			f.write(line)

	# Copy database
	os.mkdir(os.path.join(args.directory, 'database'))
	n = 0
	for f in os.listdir(dbpath):
		if f.endswith('.fasta'):
			shutil.copyfile(os.path.join(dbpath, f), os.path.join(args.directory, 'database', f))
			n += 1
	if n == 0:
		logger.error('No FASTA files in database directory. Have you selected the correct directory?')
		sys.exit(2)
	if gui is not None:
		# Only suggest to edit the config file if at least one GUI dialog has been shown
		if gui.yesno('Directory initialized',
			   'Do you want to edit the configuration file now?'):
			subprocess.call(["xdg-open", os.path.join(args.directory, PIPELINE_CONF)])
	logger.info('Directory %s initialized.', args.directory)
	logger.info('Edit %s/%s, then run "cd %s && igdiscover run" to start the analysis', args.directory, PIPELINE_CONF, args.directory)
