============
Installation
============

.. _quickstart:

Quickstart
----------

Try these instructions first. If they don’t work, continue reading below.


Linux
~~~~~

On Linux, open a terminal window and run these commands. Copy and paste
each line separately::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh

Close the terminal window and open a new one. Then run::

	conda install -c bioconda igdiscover python=3.5
	igdiscover --version

If the last command outputs the current version of IgDiscover, it worked! Skip
to the next section if you have problems.


OS X/macOS
~~~~~~~~~~

On macOS, `follow the miniconda <https://conda.io/docs/install/quick.html#os-x-miniconda-install>`_
installation instructions. Then open a terminal window and run::

	conda install -c bioconda igdiscover python=3.5
	igdiscover --version

If the last command outputs the current version of IgDiscover, it worked!
Keep reading if you have problems.


.. _simple-installation:

Conda
-----

IgDiscover is written in Python 3 and runs on Linux and OS X. Be aware that running on OS X is not
as well tested as on Linux, yet.

IgDiscover is easy to install because the program itself and its dependencies are available as
`Conda <http://conda.pydata.org/docs/>`_ packages. Conda makes it possible to install IgDiscover
with a single command, but Conda itself needs to be installed first. Since IgDiscover depends on
many programs that are not available in most Linux distributions (or OS X for that matter), we
recommend using Conda, but :ref:`non-Conda installation instructions <manual-installation>` are
also available.

First, `download and install miniconda <http://conda.pydata.org/docs/install/quick.html>`_.
That page has instructions for both Linux and OS X. If you are on Linux, either follow the
instructions there or try these commands::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh

Accept the license, and when you get the question “Do you wish the installer to
prepend the Miniconda3 install location to PATH in your /home/.../.bashrc”,
answer ``yes``. If you don’t, you need to note down the ``export PATH=...``
that the installer prints somewhere and type it in every time before using
Conda or IgDiscover.

.. note::
    If you use ``zsh`` instead of ``bash`` (applies to Bio-Linux, for example),
    then you need to manually add the ``export PATH=...`` line to the end of
    the file ``/home/your-user-name/.zshrc`` instead. The conda installer does
    not do this for you!

    Alternatively, change your default shell to bash by running
    ``chsh -s /bin/bash``.

Next, close the terminal window and open a new one. Make sure that Conda is
installed by typing in ``conda list``. If everything is working, you should see
a list of some packages andd IgDiscover can now be installed.

Install IgDiscover from `the bioconda channel <https://bioconda.github.io/bioconda>`_
with this command::

	conda install -c bioconda igdiscover python=3.5

The ``python=3.5`` is needed because Conda comes with a default of Python 3.6,
but bioconda offers only packages for Python 3.5.

Then make sure that it works by running ::

	igdiscover --version

and you should see the version number of IgDiscover.

You should now :ref:`run IgDiscover on the test data set <test>`.


If you cannot or do not want to use Conda, you need to follow the
:ref:`manual installation instructions <manual-installation>`.


Troubleshooting on OS X
-----------------------

If you get the error ::

    ValueError: unknown locale: UTF-8

Then follow `these instructions <http://conda.pydata.org/docs/troubleshooting.html#unknown-locale>`_.


Troubleshooting on Linux
------------------------

If you use conda and see an error that includes something like this::

    ImportError: .../.local/lib/python3.5/site-packages/sqt/_helpers.cpython-35m-x86_64-linux-gnu.so: undefined symbol: PyFPE_jbuf

Or you see any error that mentions a ``.local/`` directory, then a previous
installation of IgDiscover is interfering with the conda installation.

The easiest way to solve this problem is to delete the directory ``.local/`` in
your home directory, see also :ref:`how to remove IgDiscover from a Linux
system <removing-igdiscover>`.

Development version
-------------------

To install IgDiscover directly from the most recent source code,
:ref:`read the developer installation instructions <developer-install>`.
