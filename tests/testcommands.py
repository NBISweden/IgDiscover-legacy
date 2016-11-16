from igdiscover.__main__ import main
from nose.tools import assert_raises


def test_main():
	with assert_raises(SystemExit) as exc:
		main(['--version'])
	assert exc.exception.code == 0




