============
Installation
============

To install IgY-Pipe, you first need to install the programs that it requires:

* FastQC
* vsearch
* IgBLAST
* pear
* flash (this is optional)
* HDF5 development libraries (under Debian/Ubuntu, install the package
  'libhdf5-dev')

For the pipeline to be able to run them, the tools must be available in your
``$PATH``.

To then install IgY-Pipe, run ::

	python3 setup.py install --user

Then add the directory ``$HOME/.local/bin`` to your ``$PATH``.

The ``$IGDATA`` environment variable also needs to be set to the directory that
contains the ``internal_data`` and ``optional_file`` directories of IgBLAST.


Notes
-----

apt-get install python3 pigz libhdf5-dev

Install vsearch::
	wget -O ~/.local/bin/vsearch https://github.com/torognes/vsearch/releases/download/v1.1.0/vsearch-1.1.0-linux-x86_64
	chmod +x ~/.local/bin/vsearch

Install PEAR::
	wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
	tar xvf pear-0.9.6-bin-64.tar.gz
	mv pear-0.9.6-bin-64/pear-0.9.6-bin-64 ~/.local/bin/pear

Install IgBLAST::
	wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/ncbi-igblast-1.4.0-x64-linux.tar.gz
	tar xvf ncbi-igblast-1.4.0-x64-linux.tar.gz
	mv ncbi-igblast-1.4.0/bin/igblast? ~/.local/bin/
	mkdir ~/.local/igdata
	cd ~/.local/igdata
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file/

Add the following lines at the end of ``~/.bashrc``::
	export PATH=$HOME/.local/bin:$PATH
	export IGDATA=$HOME/.local/igdata
