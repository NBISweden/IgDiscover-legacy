import sys
from setuptools import setup
import versioneer


if sys.version_info < (3, 6):
	sys.stdout.write("At least Python 3.6 is required.\n")
	sys.exit(1)

with open('README.rst', encoding='utf-8') as f:
	long_description = f.read()

setup(
	name = 'igdiscover',
	version = versioneer.get_version(),
	cmdclass = versioneer.get_cmdclass(),
	author = 'Marcel Martin',
	author_email = 'marcel.martin@scilifelab.se',
	url = 'https://igdiscover.se/',
	description = 'Analyze antibody repertoires and discover new V genes',
	long_description = long_description,
	license = 'MIT',
	entry_points = {'console_scripts': ['igdiscover = igdiscover.__main__:main']},
	packages = ['igdiscover'],
	package_data = {'igdiscover': ['igdiscover.yaml', 'Snakefile', 'empty.aux']},
	python_requires = '>=3.6',
	install_requires = [
		'sqt>=0.8.0',
		'pandas>=0.23',
		'numpy',
		'matplotlib>=3.0',
		'seaborn>=0.9',
		'snakemake>=4.5',
		'cutadapt',
		'scipy>=1.1',
		'xopen>=0.3.5',
		'ruamel.yaml',
	],
	classifiers = [
		"Development Status :: 4 - Beta",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
