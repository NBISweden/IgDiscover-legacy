"""
Test the SisterMerger class
"""
import pandas as pd
from igdiscover.discover import SisterMerger, SisterInfo
from igdiscover.compose import SequenceMerger, SequenceInfo
from igdiscover.utils import UniqueNamer
from igdiscover.rename import PrefixDict
from nose.tools import raises

def test_0():
	merger = SisterMerger()
	assert list(merger) == []


def test_1():
	merger = SisterMerger()
	group = pd.DataFrame([1, 2, 3])
	info = SisterInfo('ACCGGT', False, 'name1', group)
	merger.add(info)
	assert list(merger) == [info]


def test_2():
	merger = SisterMerger()
	group = pd.DataFrame([1, 2, 3])
	merger.add(SisterInfo('ACCGGT', False, 'name1', group))
	merger.add(SisterInfo('ACCGGT', False, 'name2', group))
	sisters = list(merger)
	assert len(sisters) == 1
	assert sisters[0].sequence == 'ACCGGT'
	assert not sisters[0].requested
	assert sisters[0].name == 'name1;name2'


def test_requested():
	merger = SisterMerger()
	group = pd.DataFrame([1, 2, 3])
	merger.add(SisterInfo('ACCGGT', True, 'name1', group))
	merger.add(SisterInfo('ACCGGT', False, 'name2', group))
	sisters = list(merger)
	assert sisters[0].requested


def test_prefix():
	merger = SisterMerger()
	group = pd.DataFrame([1, 2, 3])
	merger.add(SisterInfo('ACCGGTAACGT', True, 'name1', group))
	merger.add(SisterInfo('ACCGGT', False, 'name2', group))
	info2 = SisterInfo('TGATACC', False, 'name3', group)
	merger.add(info2)
	sisters = list(merger)
	assert len(sisters) == 2
	assert sisters[0].sequence == 'ACCGGTAACGT'
	assert sisters[0].name == 'name1;name2'
	assert sisters[1] == info2


def test_with_N():
	merger = SisterMerger()
	group = pd.DataFrame([1, 2, 3])
	merger.add(SisterInfo('ACCNGTAANGT', True, 'name1', group))
	merger.add(SisterInfo('ANCGGT', False, 'name2', group))
	info2 = SisterInfo('TGATACC', False, 'name3', group)
	merger.add(info2)
	sisters = list(merger)
	assert len(sisters) == 2
	assert sisters[0].sequence == 'ACCGGTAANGT'
	assert sisters[0].name == 'name1;name2'
	assert sisters[1] == info2


def test_unique_namer():
	namer = UniqueNamer()
	assert namer('Name') == 'Name'
	assert namer('AnotherName') == 'AnotherName'
	assert namer('Name') == 'NameA'
	assert namer('Name') == 'NameB'
	assert namer('YetAnotherName') == 'YetAnotherName'
	assert namer('Name') == 'NameC'
	assert namer('NameC') == 'NameCA'


def test_sequence_merger_withCDR3():
	merger = SequenceMerger(max_differences=1)
	infos = [
		SequenceInfo('ACGTTA', 'Name1', 15),
		SequenceInfo('ACGTTAT', 'Name2', 100),  # kept because it is longer
		SequenceInfo('ACGCCAT', 'Name3', 15),  # kept because edit distance > 1
		SequenceInfo('ACGGTAT', 'Name5', 120),  # kept because more CDR3s
	]
	merger.add(infos[0]); merged = list(merger)
	assert len(merged) == 1 and merged[0] == infos[0]
	merger.add(infos[0]); merged = list(merger)
	assert len(merged) == 1 and merged[0] == infos[0]
	merger.add(infos[1]); merged = list(merger)
	assert len(merged) == 1 and merged[0] == infos[1]
	merger.add(infos[2]); merged = list(merger)
	assert len(merged) == 2 and merged[0] == infos[1] and merged[1] == infos[2]
	merger.add(infos[3]); merged = list(merger)
	assert len(merged) == 2 and merged[0] == infos[2] and merged[1] == infos[3]


def test_sequence_merger_prefix():
	merger = SequenceMerger(max_differences=1)
	infos = [
		SequenceInfo('AAATAA', 'Name1', 117),
		SequenceInfo('AAACAAG', 'Name2', 10),
	]
	merger.add(infos[0])
	merger.add(infos[1])
	merged = list(merger)
	assert len(merged) == 1
	assert merged[0] == infos[0]


def test_merger_check_all_previous():
	merger = SequenceMerger(max_differences=1)
	infos = [
		SequenceInfo('ATAAAA', 's1', 11),
		SequenceInfo('AAGAAA', 's2', 12),
		SequenceInfo('AAACAA', 's3', 13),
		SequenceInfo('AAAAAA', 'Final', 200)
	]
	for info in infos:
		merger.add(info)
	merged = list(merger)
	assert len(merged) == 1
	assert merged[0] == infos[3]



class TestPrefixDict:
	def setup(self):
		self.pd = PrefixDict([
			('AAACCT', 7),
			('AGAAA', 11),
			('AGAAC', 13)
		])

	@raises(KeyError)
	def test_ambiguous(self):
		self.pd['AGAA']

	@raises(KeyError)
	def test_missing(self):
		self.pd['TTAAG']

	def test_existing(self):
		assert self.pd['AAACCT'] == 7
		assert self.pd['AAACCT'] == 7
		assert self.pd['AAA'] == 7
		assert self.pd['AGAAAT'] == 11
		assert self.pd['AGAACT'] == 13
