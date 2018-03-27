import os
import sys
from tempfile import TemporaryDirectory
import pytest

from igdiscover.__main__ import main
from .utils import datapath, resultpath, files_equal


def run(args, expected):
	"""
	Run IgDiscover, redirecting stdout to a temporary file.
	Then compare the output with the contents of an expected file.
	"""
	with TemporaryDirectory() as td:
		outpath = os.path.join(td, 'output')
		print('Running:', ' '.join(args))
		with open(outpath, 'w') as f:
			old_stdout = sys.stdout
			sys.stdout = f
			main(args)
			sys.stdout = old_stdout
		assert files_equal(expected, outpath)


def test_main():
	with pytest.raises(SystemExit) as exc:
		main(['--version'])
	assert exc.value.code == 0


def test_group_by_barcode_only():
	args = ['group', '-b', '4', datapath('ungrouped.fasta')]
	run(args,  resultpath('grouped-by-barcode-only.fasta'))


def test_group_by_pseudo_cdr3():
	args = ['group', '-b', '4', '--pseudo-cdr3=-5:-2', '--trim-g', datapath('ungrouped.fasta')]
	run(args,  resultpath('grouped.fasta'))


def test_group_by_pseudo_cdr3_barcode_at_end():
	args = ['group', '-b', '-4', '--pseudo-cdr3=1:3', datapath('ungrouped.fasta')]
	run(args, resultpath('grouped2.fasta'))


def test_clusterplot(tmpdir):
	main(['clusterplot', '-m', '10', datapath('clusterplot.tab.gz'), str(tmpdir)])
	assert tmpdir.join('IGHV1-1801.png').check()
