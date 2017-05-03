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
       python3 -m pip install -e .

If you are a developer, you may also want to look into Conda environments and
use those.

Whenever you want to update the software::

    cd IgDiscover
    git pull

It may also be necessary to repeat the ``python3 -m pip install -e .`` step.


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

.. _removing-igdiscover:

Removing IgDiscover from a Linux system
---------------------------------------

If you have been playing around with different installation methods (``pip``,
``conda``, ``git``, ``python3 setup.py install`` etc.) you may have multiple
copies of IgDiscover on your system and you will likely run into problems
on updates. Here is a list you can follow in order to get rid of the
installations as preparation for a clean re-install. *Do not* add ``sudo`` to
the commands below if you get permission problems, unless explicitly told to do
so! If one of the steps does not work, that is fine, just continue.

1. Delete miniconda: Run the command ``which conda``. The output will be
   something like ``/home/myusername/miniconda3/bin/conda``. The part before
   ``bin/conda`` is the miniconda installation directory. Delete that folder. In
   this case, you would need to delete ``miniconda3`` in ``/home/myusername``.
2. Run ``pip3 uninstall igdiscover``. If this runs successfully and prints some
   messages about removing files, then *repeat the same command*! Do this
   until you get a message telling you that the package cannot be uninstalled
   because it is not installed.
3. Repeat the previous step, but with ``pip3 uninstall sqt``.
4. If you have a directory named ``.local`` within your home directory, you may
   want to rename it: ``mv .local dot-local-backup`` You can also delete it, but
   there is a small risk that other software (not IgDiscover) uses that
   directory. The directory is hidden, so a normal ``ls`` will not show it.
   Use ``ls -la`` while in your home directory to see it.
5. If you have ever used ``sudo`` to install IgDiscover, you may have an
   installation in ``/usr/local/``. You can try to remove it with
   ``sudo pip3 uninstall igdiscover``.
6. Delete the cloned Git repository if you have one. This is the directory in
   which you run ``git pull``.

Finally, you can follow the normal installation instructions and then the
developer installation instructions.
