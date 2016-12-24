import sys
from setuptools import setup
import versioneer


if sys.version_info < (3, 4):
	sys.stdout.write("At least Python 3.4 is required.\n")
	sys.exit(1)

with open('README.rst') as f:
	long_description = f.read()

setup(
	name = 'igdiscover',
	version = versioneer.get_version(),
	cmdclass = versioneer.get_cmdclass(),
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = 'https://igdiscover.readthedocs.io/',
	description = 'Analyze antibody repertoires and discover new V genes',
	long_description = long_description,
	license = 'MIT',
	entry_points = {'console_scripts': ['igdiscover = igdiscover.__main__:main']},
	packages = ['igdiscover'],
	package_data = {'igdiscover': ['igdiscover.yaml', 'Snakefile', 'empty.aux']},
	install_requires = [
		'sqt>=0.8.0',
		'pandas>=0.16.2',
		'numpy',
		'matplotlib>=1.5.0',
		'snakemake>=3.9.0',
		'cutadapt',
		'seaborn>=0.6.0',
		'scipy>=0.16.1',
		'xopen>=0.1.1',
		'ruamel.yaml',
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
