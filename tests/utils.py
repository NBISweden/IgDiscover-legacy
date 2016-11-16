import sys
import os
import subprocess
from contextlib import contextmanager
from io import StringIO


@contextmanager
def capture_stdout():
	sio = StringIO()
	old_stdout = sys.stdout
	sys.stdout = sio
	yield sio
	sys.stdout = old_stdout


def datapath(path):
	return os.path.join(os.path.dirname(__file__), 'data', path)


def resultpath(path):
	return os.path.join(os.path.dirname(__file__), 'results', path)


def files_equal(path1, path2):
	return subprocess.call(['diff', '-u', path1, path2]) == 0

