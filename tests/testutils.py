from io import StringIO
import pkg_resources
from igdiscover.utils import has_stop
from igdiscover.utils import Config


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


def test_config():
	empty_config = Config(file=StringIO('{}'))
	packaged_config = Config(file=pkg_resources.resource_stream('igdiscover', Config.DEFAULT_PATH))
	# force library name to be equal since it is dynamically determined
	empty_config.library_name = packaged_config.library_name = 'nolib'
	e = empty_config.__dict__
	p = packaged_config.__dict__
	assert e.keys() == p.keys()
	for k in e:
		assert e[k] == p[k], '{}: {} vs {}'.format(k, e[k], p[k])
	# assert empty_config == packaged_config
