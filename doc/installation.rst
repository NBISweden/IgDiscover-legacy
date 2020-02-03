============
Installation
============

IgDiscover is written in Python 3 and is developed on Linux. The tool also
runs on macOS, but is not as well tested on that platform.

For installation on either system, we recommend that you follow the instructions
below, which will first explain how to install the `Conda <http://conda.io/docs/>`_
package manager. IgDiscover is available as a
Conda-package from `the bioconda channel <https://bioconda.github.io/bioconda>`_.
Using Conda will make the installation easy because all dependencies are also
available as Conda packages and can thus be installed automatically along with
IgDiscover.

There are also :ref:`non-Conda installation instructions <manual-installation>`
if you cannot use Conda.


.. _install-with-conda:

Installing IgDiscover with Conda
--------------------------------

1. Install `Conda`_ by following the `conda installation
   instructions <https://conda.io/docs/user-guide/install/>`_
   as appropriate for your system. You will need to choose between a “Miniconda”
   and “Anaconda” installation. We recommend Miniconda as the download is
   smaller. If you are in a hurry, these two commands are usually sufficient to
   install Miniconda on Linux (read the linked document for macOS instructions)::

       wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
       bash Miniconda3-latest-Linux-x86_64.sh

   When the installer asks you about modifying the ``PATH`` in your ``.bashrc``
   file, answer ``yes``.

2. Close the terminal window and open a new one. Then test whether conda is
   installed correctly by running ::

       conda --version

   If you see the conda version number, it worked.
3. Set up Conda so that it can access the
   `bioconda channel <https://bioconda.github.io/>`_.
   For that, follow `the instructions on the bioconda
   website <https://bioconda.github.io/#set-up-channels>`_
   or simply run these commands::

       conda config --add channels defaults
       conda config --add channels bioconda
       conda config --add channels conda-forge
4. Install IgDiscover with this command::

       conda create -n igdiscover igdiscover

   This will create a new so-called “environment” for IgDiscover (retry if it fails). **Whenever you
   want to run IgDiscover, you will need to activate the environment with this
   command**::

       source activate igdiscover

5. Make sure you have activated the ``igdiscover`` environment.
   Then test whether IgDiscover is correctly installed with this command::

       igdiscover --version

   If you see the version number of IgDiscover, it worked! If an error message appears that says 
   "The 'networkx' distribution was not found and is required by snakemake", install networkx manually with::
      
       pip install networkx==2.1

   Then retry to check the igdiscover version.

6. You can now :ref:`run IgDiscover on the test data set <test>` to familiarize
   yourself with how it works.


.. _troubleshooting:


Troubleshooting on Linux
------------------------

If you use ``zsh`` instead of ``bash`` (applies to Bio-Linux, for example),
the ``$PATH`` environment variable will not be setup correctly by the
Conda installer. The miniconda installer adds a line ``export PATH=...`` to the
to the end of your ``/home/your-user-name/.bashrc`` file. Copy that line from
the file and add it to the end of the file ``/home/your-user-name/.zshrc``
instead.

Alternatively, change your default shell to bash by running
``chsh -s /bin/bash``.

If you use conda and see an error that includes something like this::

    ImportError: .../.local/lib/python3.5/site-packages/sqt/_helpers.cpython-35m-x86_64-linux-gnu.so: undefined symbol: PyFPE_jbuf

Or you see any error that mentions a ``.local/`` directory, then a previous
installation of IgDiscover is interfering with the conda installation.

The easiest way to solve this problem is to delete the directory ``.local/`` in
your home directory, see also :ref:`how to remove IgDiscover from a Linux
system <removing-igdiscover>`.


Troubleshooting on macOS
------------------------

If you get the error ::

    ValueError: unknown locale: UTF-8

Then follow `these instructions <https://conda.io/docs/user-guide/troubleshooting.html#macos-error-valueerror-unknown-locale-utf-8>`_.


Development version
-------------------

To install IgDiscover directly from the most recent source code,
:ref:`read the developer installation instructions <developer-install>`.
