from io import StringIO
import pkg_resources
import pytest

from igdiscover.utils import (has_stop, validate_fasta, FastaValidationError, find_overlap,
    merge_overlapping, UniqueNamer)
from igdiscover.cli.config import Config


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
    with pkg_resources.resource_stream('igdiscover', Config.DEFAULT_PATH) as file:
        packaged_config = Config(file=file)
    # force library name to be equal since it is dynamically determined
    empty_config.library_name = packaged_config.library_name = 'nolib'
    e = empty_config.__dict__
    p = packaged_config.__dict__
    assert_dicts_equal(e, p)
    # assert empty_config == packaged_config


def test_validate_empty_record():
    with pytest.raises(FastaValidationError):
        validate_fasta('tests/data/empty-record.fasta')


def test_validate_duplicate_name():
    with pytest.raises(FastaValidationError):
        validate_fasta('tests/data/duplicate-name.fasta')


def test_validate_duplicate_sequence():
    with pytest.raises(FastaValidationError):
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


def test_unique_namer():
    namer = UniqueNamer()
    assert "ABC_S1234" == namer("ABC_S1234")
    assert "ABC_S1234A" == namer("ABC_S1234")
    assert "X" == namer("X")
    assert "ABC_S1234B" == namer("ABC_S1234")
