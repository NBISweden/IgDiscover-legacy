.. _develop:

Development
===========


* `Source code <https://bitbucket.org/igdiscover/igdiscover/>`_
* `Report an issue <https://bitbucket.org/igdiscover/igdiscover/issues?status=new&status=open>`_


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


