from io import StringIO
import pkg_resources
from nose.tools import raises
from igdiscover.utils import has_stop, validate_fasta, FastaValidationError, find_overlap, merge_overlapping
from igdiscover.config import Config


def test_has_stop():
	assert has_stop('TAA')
	assert has_stop('TAG')
	assert has_stop('TGA')
	assert has_stop('GGGTGA')

	assert not has_stop('GGG')
	assert not has_stop('TAT')
	assert not has_stop('GGGT')
	assert not has_stop('GGGTA')
	assert not has_stop('TATTG')


def assert_dicts_equal(expected, actual):
	assert expected.keys() == actual.keys()
	for k in expected:
		if hasattr(expected[k], 'keys'):
			assert hasattr(actual[k], 'keys')
			assert_dicts_equal(expected[k], actual[k])
		else:
			assert expected[k] == actual[k], '{}: {} vs {}'.format(k, expected[k], actual[k])


def test_config():
	empty_config = Config(file=StringIO('{}'))
	packaged_config = Config(file=pkg_resources.resource_stream('igdiscover', Config.DEFAULT_PATH))
	# force library name to be equal since it is dynamically determined
	empty_config.library_name = packaged_config.library_name = 'nolib'
	e = empty_config.__dict__
	p = packaged_config.__dict__
	assert_dicts_equal(e, p)
	# assert empty_config == packaged_config


@raises(FastaValidationError)
def test_validate_empty_record():
	validate_fasta('tests/data/empty-record.fasta')


@raises(FastaValidationError)
def test_validate_duplicate_name():
	validate_fasta('tests/data/duplicate-name.fasta')


@raises(FastaValidationError)
def test_validate_duplicate_sequence():
	validate_fasta('tests/data/duplicate-sequence.fasta')


def test_find_overlap():
	assert find_overlap('', '') is None
	assert find_overlap('A', '') is None
	assert find_overlap('ABC', 'X') is None
	assert find_overlap('X', 'ABC') is None
	assert find_overlap('A', 'A') == 0
	assert find_overlap('ABCD', 'A') == 0
	assert find_overlap('A', 'ABC') == 0
	assert find_overlap('AB', 'BD') == 1
	assert find_overlap('ABCDE', 'CDE') == 2
	assert find_overlap('ABCDEFGH', 'CDE') == 2
	assert find_overlap('CDE', 'XABCDEFG') == -3
	assert find_overlap('EFGHI', 'ABCDEFG') == -4


def test_merge_overlapping():
	assert merge_overlapping('', '') is None  # TODO
	assert merge_overlapping('ABC', 'DEF') is None
	assert merge_overlapping('HELLOW', 'LOWORLD') == 'HELLOWORLD'
	assert merge_overlapping('LOWORLD', 'HELLOW') == 'HELLOWORLD'
	assert merge_overlapping('HELLOWORLD', 'LOWO') == 'HELLOWORLD'
	assert merge_overlapping('LOWO', 'HELLOWORLD') == 'HELLOWORLD'
