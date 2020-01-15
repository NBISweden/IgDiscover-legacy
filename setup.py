import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 6):
    sys.stdout.write("At least Python 3.6 is required.\n")
    sys.exit(1)

with open("README.rst", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="igdiscover",
    setup_requires=["setuptools_scm"],  # Support pip versions that don't know about pyproject.toml
    use_scm_version=True,
    author="Marcel Martin",
    author_email="marcel.martin@scilifelab.se",
    url="https://igdiscover.se/",
    description="Analyze antibody repertoires and discover new V genes",
    long_description=long_description,
    license="MIT",
    entry_points={"console_scripts": ["igdiscover = igdiscover.__main__:main"]},
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={"igdiscover": ["igdiscover.yaml", "Snakefile", "empty.aux"]},
    python_requires=">=3.6",
    install_requires=[
        "pandas>=0.23",
        "numpy",
        "matplotlib>=3.0",
        "seaborn>=0.9",
        "snakemake>=4.5",
        "cutadapt>=2.5",
        "scipy>=1.1",
        "xopen>=0.8.2",
        "ruamel.yaml",
        "dnaio",
        "tinyalign>=0.2",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
