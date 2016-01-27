============
Installation
============

IgY-Pipe is written in Python 3 and needs a Linux distribution to run. It may
even work on OS X and even Windows, but this has not been tested.

IgY-Pipe depends on many programs that are not available in most Linux
distributions. To simplify the installation process, we provide everything as
`conda <http://conda.pydata.org/docs/>`_ packages. Conda makes it possible to
install IgY-Pipe with a single command, but conda itself needs to be installed
first. If you do not want to use conda, follow the :ref:`manual installation
instructions <manual-installation>`.


Simple installation with conda
------------------------------

First, `download and install miniconda <http://conda.pydata.org/docs/install/quick.html#linux-miniconda-install>`_.
Either follow the instructions on the linked page or try these commands::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh

Then close the terminal window and open a new one. Make sure that conda is
installed by running ``conda list``. If everything is working, you should see
a list of some packages. If it works, installing IgY-Pipe is now just a matter
of typing in a single command.

Try to install IgY-Pipe from the bioconda channel with this command::

	conda install -c bioconda igypipe

If this does not work, then we have not yet made IgY-Pipe available online.
In this case, you need to get the conda package as a ``.tar.gz`` file from us.
Then install it with this command (replace the file name with the actual file name)::

	conda install -c bioconda igypipe-0.x.y-0.tar.gz

Make sure that it works by running ::

	igypipe --version

and you should see the version number of IgY-Pipe.


.. _manual-installation:

Manual installation without conda
---------------------------------

Non-Python dependencies need to be installed by hand. When IgY-Pipe itself is
installed in the last step, the Python dependencies will automatically be
installed.

Install dependencies
~~~~~~~~~~~~~~~~~~~~

The dependencies are: MUSCLE, FastQC, VSEARCH, IgBLAST, PEAR, and --
optionally -- flash.

1. Install Python 3. It most likely is already installed on your system, but
   in Debian/Ubuntu, you can get it with ::

	sudo apt-get install python3

2. Create the directory where binaries will be installed. We assume
   ``$HOME/.local/bin`` here, but this can be anywhere as long as they are in
   your ``$PATH``. ::

	mkdir -p ~/.local/bin

   Add this line to the end of your ``~/.bashrc`` file::

	export PATH=$HOME/.local/bin:$PATH

   Then either start a new shell or run ``source ~/.bashrc`` to get the changes.

3. Install MUSCLE. This is available as a package in Ubuntu::

	sudo apt-get install muscle

   If your distribution does not have a 'muscle' package or if are not allowed
   to run ``sudo``::

	wget -O - http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz | tar xz
	mv muscle3.8.31_i86linux64 ~/.local/bin/

4. Install VSEARCH::

	wget -O ~/.local/bin/vsearch https://github.com/torognes/vsearch/releases/download/v1.1.0/vsearch-1.1.0-linux-x86_64
	chmod +x ~/.local/bin/vsearch

5. Install PEAR::

	wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
	tar xvf pear-0.9.6-bin-64.tar.gz
	mv pear-0.9.6-bin-64/pear-0.9.6-bin-64 ~/.local/bin/pear

6. Install IgBLAST::

	wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.4.0/ncbi-igblast-1.4.0-x64-linux.tar.gz
	tar xvf ncbi-igblast-1.4.0-x64-linux.tar.gz
	mv ncbi-igblast-1.4.0/bin/igblast? ~/.local/bin/

   IgBLAST requires some data files that must be downloaded separately. The
   following commands put the files into ``~/.local/igdata``::

	mkdir ~/.local/igdata
	cd ~/.local/igdata
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/
	wget -r -nH --cut-dirs=4 ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file/

   Also, you must set the ``$IGDATA`` environment variable to point to the
   directory with data files. Add this line to your ``~/.bashrc``::

	export IGDATA=$HOME/.local/igdata

   Then run ``source ~/.bashrc`` to get the changes.

7. Install FastQC. Under Debian/Ubuntu::

	sudo apt-get install fastqc

   If you install it by hand, make sure you end up with a ``fastqc`` binary in
   the ``$PATH``.

8. (Optionally, install flash. Not described here.)


Install IgY-Pipe
~~~~~~~~~~~~~~~~

If you have an ``igypipe-....tar.gz`` file, then unpack it, change into the
created directory and run ::

	python3 setup.py install --user

If you do not have a ``.tar.gz`` file, install IgY-Pipe with the Python
package manager ``pip``, which will download IgY-Pipe before installing it::

	pip3 install --user igypipe

Both commands also install all remaining dependencies. The ``--user`` option
instructs both commands to install everything into ``$HOME/.local``.

Finally, check the installation with ::

	igypipe --version

and you should see the version number of IgY-Pipe.
