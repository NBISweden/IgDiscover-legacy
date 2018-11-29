.. _manual-installation:

Manual installation
===================

IgDiscover requires quite a few other software tools that are not included in most Linux
distributions (or mac OS) and which are also not available from the Python packaging
index (PyPI) because they are not Python tools. If you do not use the :ref:`recommended simple
installation instructions via Conda <install-with-conda>`, you need to install those non-Python
dependencies manually. Regular Python dependencies are automatically pulled in when IgDiscover
itself is installed in the last step with the ``pip install`` command. The instructions below are
written for Linux and require modifications if you want to try this on OS X.

.. note::
    We recommend the much simpler :ref:`installation via Conda <install-with-conda>`
    instead of using the instructions in this section.


Install non-Python dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dependencies are: MUSCLE, IgBLAST, PEAR, and -- optionally -- flash.

1. Install Python 3.6 or newer. It may already be installed on your system. On
   a sufficiently recent Debian or Ubuntu, you can get it with ::

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

   If your distribution does not have a 'muscle' package or if you are not allowed
   to run ``sudo``::

	wget -O - http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz | tar xz
	mv muscle3.8.31_i86linux64 ~/.local/bin/

4. Install PEAR::

	wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
	tar xvf pear-0.9.6-bin-64.tar.gz
	mv pear-0.9.6-bin-64/pear-0.9.6-bin-64 ~/.local/bin/pear

5. Install IgBLAST::

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

7. Optionally, install flash::

	wget -O FLASH-1.2.11.tar.gz http://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz/download
	tar xf FLASH-1.2.11.tar.gz
	cd FLASH-1.2.11
	make
	mv flash ~/.local/bin/


Install IgDiscover
~~~~~~~~~~~~~~~~~~

Install IgDiscover with the Python package manager ``pip``, which will download and install IgDiscover and its
dependencies::

	pip3 install --user igdiscover

Both commands also install all remaining dependencies. The ``--user`` option
instructs both commands to install everything into ``$HOME/.local``.

Finally, check the installation with ::

	igdiscover --version

and you should see the version number of IgDiscover.

You should now :ref:`run IgDiscover on the test data set <test>`.
