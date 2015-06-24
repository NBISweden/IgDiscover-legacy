import sys
from setuptools import setup
from igypipe import __version__

if sys.version_info < (3, 3):
	sys.stdout.write("At least Python 3.3 is required.\n")
	sys.exit(1)

setup(
	name = 'igypipe',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = '',
	description = '',
	#license = 'MIT',
	packages = ['igypipe'],
	scripts = ['bin/igypipe'],
	install_requires = [
		'sqt>=0.4.4',
		'pandas',
		'numpy',
		'matplotlib',
		'snakemake',
		'cutadapt'
	],
	classifiers = [
		"Development Status :: 2 - Pre-Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
#		"License :: OSI Approved :: MIT License", (license undecided)
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
