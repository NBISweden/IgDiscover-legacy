import sys
from setuptools import setup

# set __version__
exec(next(open('igypipe/__init__.py')))

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
	package_data = {'igypipe': ['pipeline.conf', 'Snakefile', 'empty.aux',
		'databases/*/*.fasta']},
	install_requires = [
		'sqt>=0.4.6',
		'pandas>=0.16.2',
		'numpy',
		'matplotlib',
		'snakemake',
		'cutadapt',
		'seaborn>=0.6.0',
		'tables>=3.2.1'
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
