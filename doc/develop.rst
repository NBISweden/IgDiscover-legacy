.. _develop:

Development
===========


* `Source code <https://github.com/NBISweden/IgDiscover/>`_
* `Report an issue <https://github.com/NBISweden/IgDiscover/issues>`_


.. _developer-install:

Installing the development version
----------------------------------

To use the most recent IgDiscover version from Git, first :ref:`follow the
regular installation instructions <quickstart>`. Then do the following:

1. Uninstall the igdiscover package with::

       conda uninstall igdiscover

   This removes igdiscover itself, but the dependencies remain (IgBLAST etc.),
   which is what we need.

2. Clone the repository::

       git clone https://github.com/NBISweden/IgDiscover.git

   (Use the git@ URL instead if you are a developer.)

3. Install IgDiscover::

       cd IgDiscover
       pip3 install -e .

If you are a developer, you may also want to look into Conda environments and
use those.

Whenever you want to update the software::

    cd IgDiscover
    git pull

It may also be necessary to repeat the ``pip install -e .`` step.


Building the documentation
--------------------------

Go to the ``doc/`` directory in the repository, then run::

    make

to build the documentation locally. Open ``_build/html/index.html`` in
a browser. The layout is different from the `version shown on
Read the Docs <https://igdiscover.readthedocs.io/>`_, but allows you to
preview any changes you may have made.


Making a release
----------------

We use `versioneer <https://github.com/warner/python-versioneer>`_ to
manage version numbers. It extracts the version number from the
most recent tag in Git. Thus, to increment the version number, create
a Git tag::

    git tag v0.5

The ``v`` prefix is mandatory.


