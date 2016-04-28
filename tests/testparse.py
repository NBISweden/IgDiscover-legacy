"""
"""
from igdiscover.parse import parse_header
from nose.tools import raises


def test():
	assert parse_header('') == ('', None, None)
	assert parse_header('abc') == ('abc', None, None)
	assert parse_header('abc ;barcode=ACG') == ('abc', None, 'ACG')
	assert parse_header('abc size=17;barcode=ACG') == ('abc', 17, 'ACG')
	assert parse_header('abc \t size=17;barcode=ACG') == ('abc', 17, 'ACG')
	assert parse_header('abc ;size=17;barcode=ACG;soze=19;baarcode=GGG') == ('abc', 17, 'ACG')
