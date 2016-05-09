import sys
from setuptools import setup

# set __version__
exec(next(open('igdiscover/__init__.py')))

if sys.version_info < (3, 3):
	sys.stdout.write("At least Python 3.3 is required.\n")
	sys.exit(1)

setup(
	name = 'igdiscover',
	version = __version__,
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = '',
	description = '',
	license = 'MIT',
	entry_points = {'console_scripts': ['igdiscover = igdiscover.__main__:main']},
	packages = ['igdiscover'],
	package_data = {'igdiscover': ['igdiscover.yaml', 'Snakefile', 'empty.aux']},
	install_requires = [
		'pysam==0.8.4',
		'sqt>=0.6.1',
		'pandas>=0.16.2',
		'numpy',
		'matplotlib>=1.5.0',
		'snakemake',
		'cutadapt',
		'seaborn>=0.6.0',
		'scipy>=0.16.1',
		'PyYAML',
	],
	classifiers = [
		"Development Status :: 3 - Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
