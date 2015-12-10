"""
Test the SisterMerger class
"""
import pandas as pd
from igypipe.discover import SisterMerger, SisterInfo
from igypipe.compose import Merger, SequenceInfo

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
