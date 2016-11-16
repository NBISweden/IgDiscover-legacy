from igdiscover.__main__ import main
from nose.tools import assert_raises
from .utils import datapath, resultpath, capture_stdout, files_equal
from tempfile import TemporaryDirectory
import os
import sys

def test_main():
	with assert_raises(SystemExit) as exc:
		main(['--version'])
	assert exc.exception.code == 0


def test_group():
	with TemporaryDirectory() as td:
		path = os.path.join(td, 'grouped.fasta')
		cmd = ['group', '-b', '4', '--pseudo-cdr3-range=-5:-2', '--trim-g', datapath('ungrouped.fasta')]
		print('Running:', ' '.join(cmd))
		with open(path, 'w') as f:
			old_stdout = sys.stdout
			sys.stdout = f
			main(cmd)
			sys.stdout = old_stdout
		assert files_equal(resultpath('grouped.fasta'), path)
