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
import dnaio
from xopen import xopen

from . import CommandLineError
from ..config import Config

try:
    import tkinter as tk
    from tkinter import messagebox, filedialog, TclError
except ImportError:
    tk = None


logger = logging.getLogger(__name__)

do_not_show_cpustats = 1


class GuiCancelledError(CommandLineError):
    pass


def add_arguments(parser):
    parser.add_argument('--database', '--db', metavar='PATH', default=None,
        help='Directory with V.fasta, D.fasta and J.fasta files. If not given, a dialog is shown.')
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
        raise CommandLineError(f'Could not open {path!r}: {e}')


def read_and_repair_fasta(path):
    """
    Read a FASTA file and make sure it is suitable for use with makeblastdb.
    It repairs the following issues:

    - If a record is empty, it is skipped
    - If a record name occurs more than once, the second record name gets a suffix
    - If a sequence occurs more than once, occurrences after the first are skipped
    """
    with dnaio.open(path) as sr:
        records = list(sr)

    names = set()
    sequences = dict()
    for r in records:
        r.sequence = r.sequence.upper()
        if len(r.sequence) == 0:
            logger.info("Record %r is empty, skipping it.", r.name)
            continue
        name = r.name.split()[0]
        i = 0
        while name in names:
            i += 1
            name = r.name + '_{}'.format(i)
        if name != r.name.split()[0]:
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
    run_init(**vars(args))


def run_init(
    directory,
    database: str,
    reads1=None,
    single_reads=None,
):
    if ' ' in directory:
        raise CommandLineError('The name of the analysis directory must not contain spaces')

    if os.path.exists(directory):
        raise CommandLineError(f'The target directory {directory!r} already exists.')

    # If reads files or database were not given, initialize the GUI
    if (reads1 is None and single_reads is None) or database is None:
        try:
            gui = TkinterGui()
        except (ImportError, TclError) as e:
            raise CommandLineError(f'GUI cannot be started: {e}\n'
                'Please provide reads1 file and database directory on command line.'
            )
    else:
        gui = None

    # Find out whether data is paired or single
    assert not (reads1 and single_reads)
    if reads1 is None and single_reads is None:
        paired = gui.yesno('Paired end or single-end reads',
            'Are your reads paired and need to be merged?\n\n'
            'If you answer "Yes", next select the FASTQ files '
            'with the first of your paired-end reads.\n'
            'If you answer "No", next select the FASTA or FASTQ '
            'file with single-end reads.')
        if paired is None:
            raise GuiCancelledError()
    else:
        paired = bool(reads1)

    # Assign reads1 and (if paired) also reads2
    if paired:
        if reads1 is not None:
            try_open(reads1)
        else:
            reads1 = gui.reads1_path()
            if not reads1:
                raise GuiCancelledError()
        reads2 = guess_paired_path(reads1)
        if reads2 is None:
            raise CommandLineError('Could not determine second file of paired-end reads')
    else:
        if single_reads is not None:
            reads1 = single_reads
            try_open(reads1)
        else:
            reads1 = gui.single_reads_path()
            if not reads1:
                raise GuiCancelledError()

    if database is not None:
        dbpath = database
    else:
        # TODO as soon as we distribute our own database files, we can use this:
        # database_path = pkg_resources.resource_filename('igdiscover', 'databases')
        databases_path = None
        dbpath = gui.database_path(databases_path)
        if not dbpath:
            raise GuiCancelledError()

    database = dict()
    for g in ['V', 'D', 'J']:
        path = os.path.join(dbpath, g + '.fasta')
        if not os.path.exists(path):
            raise CommandLineError(
                f'The database directory {dbpath!r} must contain the three files '
                'V.fasta, D.fasta and J.fasta. A dummy D.fasta is necessary even '
                'if analyzing light chains (see manual)'
            )
        database[g] = list(read_and_repair_fasta(path))

    # Create the directory
    try:
        os.mkdir(directory)
    except OSError as e:
        raise CommandLineError(e)

    def create_symlink(readspath, dirname, target):
        gz = '.gz' if readspath.endswith('.gz') else ''
        if not os.path.isabs(readspath):
            src = os.path.relpath(readspath, dirname)
        else:
            src = readspath
        os.symlink(src, os.path.join(dirname, target + gz))

    if paired:
        create_symlink(reads1, directory, 'reads.1.fastq')
        create_symlink(reads2, directory, 'reads.2.fastq')
    else:
        try:
            target = 'reads.' + file_type(reads1)
        except UnknownFileFormatError:
            raise CommandLineError('Cannot determine whether reads file is FASTA or FASTQ')
        create_symlink(reads1, directory, target)

    # Write the configuration file
    configuration = pkg_resources.resource_string('igdiscover', Config.DEFAULT_PATH).decode()
    with open(os.path.join(directory, Config.DEFAULT_PATH), 'w') as f:
        f.write(configuration)

    # Create database files
    database_dir = os.path.join(directory, 'database')
    os.mkdir(database_dir)
    for gene in ['V', 'D', 'J']:
        with open(os.path.join(database_dir, gene + '.fasta'), 'w') as db_file:
            for record in database[gene]:
                print('>{}\n{}'.format(record.name, record.sequence), file=db_file)

    if gui is not None:
        # Only suggest to edit the config file if at least one GUI dialog has been shown
        if gui.yesno('Directory initialized',
                'Do you want to edit the configuration file now?'):
            launch(os.path.join(directory, Config.DEFAULT_PATH))
    logger.info('Directory %s initialized.', directory)
    logger.info('Edit %s/%s, then run "cd %s && igdiscover run" to start the analysis',
        directory, Config.DEFAULT_PATH, directory)
