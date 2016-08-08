from igdiscover.utils import has_stop


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
