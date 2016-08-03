============
Installation
============

IgDiscover is written in Python 3 and needs a Linux distribution to run. It may
even work on OS X, but this has not been tested.

IgDiscover is very easy to install because IgDiscover itself and its dependencies are available as
`Conda <http://conda.pydata.org/docs/>`_ packages. Conda makes it possible to install IgDiscover with a single command,
but Conda itself needs to be installed first. Since IgDiscover depends on many programs that are not available in most
Linux distributions, we recommend using Conda, but :ref:`non-Conda installation instructions <manual-installation>` are
also available.


.. _simple-installation:

Simple installation with Conda
------------------------------

First, `download and install miniconda <http://conda.pydata.org/docs/install/quick.html#linux-miniconda-install>`_.
Either follow the instructions on the linked page or try these commands::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh

Accept the license, and when you get the question “Do you wish the installer to
prepend the Miniconda3 install location to PATH in your /home/.../.bashrc”,
answer ``yes``. If you don’t, you need to type in the line ``export PATH=...``
that the installer prints every time before using Conda or IgDiscover.

Next, close the terminal window and open a new one. Make sure that Conda is
installed by typing in ``conda list``. If everything is working, you should see
a list of some packages. If it works, IgDiscover can now be installed.

Install IgDiscover from the bioconda channel with this command::

	conda install -c bioconda -c marcelm igdiscover

Then make sure that it works by running ::

	igdiscover --version

and you should see the version number of IgDiscover.

You should now :ref:`run IgDiscover on the test data set <test>`.


If you cannot or do not want to use Conda, you need to follow the
:ref:`manual installation instructions <manual-installation>`.
