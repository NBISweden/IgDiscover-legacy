"""
Create and initialize a new pipeline directory.
"""
import glob
import logging
import os
import os.path
import sys
import shutil
import subprocess
import pkg_resources

logger = logging.getLogger(__name__)


PIPELINE_CONF = 'pipeline.conf'


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('init', help=__doc__)
	subparser.set_defaults(func=init_command)
	subparser.add_argument('--database', '--db', metavar='PATH', default=None,
		help='Directory with IgBLAST database files. If not given, a dialog is shown.')
	subparser.add_argument('directory', help='New pipeline directory to create')
	subparser.add_argument('--reads1', default=None,
		help='File with paired-end reads (first file only). If not given, a dialog is shown.')
	return subparser


def tkinter_reads_path(directory=False):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	path = filedialog.askopenfilename(title="Choose first reads file",
		filetypes=[
			("Reads", "*.fasta *.fastq *.fastq.gz *.fasta.gz"),
			("Any file", "*")])
	# messagebox.showinfo('Chosen file', repr(path))
	return path


def tkinter_database_path(initialdir):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	path = filedialog.askdirectory(title="Choose IgBLAST database directory", mustexist=True, initialdir=initialdir)
	# messagebox.showinfo('Chosen file', repr(path))
	return path


def yesno(title, question):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	return messagebox.askyesno(title, question)


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


def guess_paired_end_second_file(path):
	"""
	Given the path to a file that contains the sequences for the first read in a
	pair, return the file that contains the sequences for the second read in a
	pair. Both files must have identical names, except that the first must have
	a '1' in its name, and the second must have a '2' at the same position.

	Return None if no second file was found or if there are too many candidates.

	>>> find_paired_end_second('file.1.fastq.gz')
	file.2.fastq.gz  # if that file exists

	TODO
	Will currently not work if files have multiple positions at which one has a
	'1' and the other a '2'.
	"""
	paths = set(glob.glob(path.replace('1', '?')))
	paths.remove(path)
	if len(paths) != 1:
		return None
	return paths.pop()


def init_command(args):
	gui = False
	if args.reads1 is not None:
		reads1 = args.reads1
	else:
		gui = True
		reads1 = tkinter_reads_path()
	if reads1 == '':
		logger.error('Cancelled')
		sys.exit(2)
	reads2 = guess_paired_end_second_file(reads1)
	if reads2 is None:
		logger.error('Could not determine second file of paired-end reads')
		sys.exit(1)

	if args.database is not None:
		dbpath = args.database
	else:
		gui = True
		databases_path = pkg_resources.resource_filename('igypipe', 'databases')
		dbpath = tkinter_database_path(databases_path)
		if dbpath == '':
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

	create_symlink(reads1, args.directory, 'reads.1.fastq')
	create_symlink(reads2, args.directory, 'reads.2.fastq')

	snakepath = pkg_resources.resource_filename('igypipe', 'Snakefile')
	os.symlink(os.path.relpath(snakepath, args.directory), os.path.join(args.directory, 'Snakefile'))

	# Write the pipeline configuration
	configuration = pkg_resources.resource_string('igypipe', PIPELINE_CONF)
	with open(os.path.join(args.directory, PIPELINE_CONF), 'wb') as f:
		f.write(configuration)

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
	if gui:
		# Only suggest to edit the config file if at least one GUI dialog has been shown
		if yesno('Directory initialized', 'Do you want to edit the configuration file now?'):
			subprocess.call(["xdg-open", os.path.join(args.directory, PIPELINE_CONF)])
	else:
		logger.info('Directory %s initialized.', args.directory)
	logger.info('Edit %s/pipeline.conf, then run "cd %s && snakemake -j" to start the pipeline', args.directory, args.directory)
