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
	entry_points = {'console_scripts': ['igypipe = igypipe.__main__:main']},
	packages = ['igypipe'],
	package_data = {'igypipe': ['pipeline.conf', 'Snakefile', 'empty.aux',
		'databases/*/*.fasta']},
	setup_requires = ['numpy'],
	install_requires = [
		'sqt>=0.6.0',
		'pandas>=0.16.2',
		'numpy',
		'matplotlib>=1.5.0',
		'snakemake',
		'cutadapt',
		'seaborn>=0.6.0',
#		'tables>=3.2.1',
		'scipy>=0.16.1',
		'PyYAML',
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
