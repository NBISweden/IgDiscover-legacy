.. _advanced:

Advanced topics
===============

IgDiscover itself does not (yet) come with all imaginable analysis facilities built into it.
However, it creates many files (mostly with tables) that can be used for custom analysis.
For example, all ``.tab`` and ``.tsv`` files (in particular ``assigned.tsv.gz`` and
``candidates.tab``)
can be opened and inspected in a spreadsheet application such as LibreOffice. From there,
you can do basic tasks such as sorting from the menu of that application.

Often, these facilities are not enough, however, and some basic understanding of the
command-line is helpful. Clearly, this is not as convenient as working in a graphical
user interface (GUI), but we do not currently have the resources to provide one for
IgDiscover. To alleviate this somewhat, we provide here instructions for a few things
that you may want to do with the IgDiscover result files.


Extract all sequences that match any database gene exactly
----------------------------------------------------------

The ``candidates.tab`` file tells you for each discovered sequence how often an *exact match*
of that sequence was found in your input reads. A high number of exact matches is a good
indication that the candidate is actually a new gene or allele. In order to find the original
reads that correspond to those matches, you can run this command in the analysis directory,
replacing ``iteration-01`` with the directory in which the ``filtered.tsv.gz`` file is
located::

    igdiscover run iteration-01/exact.tab

This command will extract all rows from ``iteration-01/filtered.tsv.gz`` for which the
``V_errors`` column is zero.


Extra configuration settings
----------------------------

Some configuration settings are not documented in the default ``igdiscover.yaml`` file
since they rarely need to be changed.

::

    # Leave empty or choose a species name supported by IgBLAST:
    # human, mouse, rabbit, rat, rhesus_monkey
    # This setting is not used anywhere except that it is passed
    # to IgBLAST using the -organism option. Since we provide IgBLAST
    # with our own gene databases, it seems this has no effect.
    species:

::

    # Which program to use for computing multiple alignments. This is used for
    # computing consens sequences.
    # Choose 'mafft', 'clustalo', 'muscle' or 'muscle-fast'.
    # 'muscle-fast' runs muscle with parameters "-maxiters 1 -diags".
    #
    #multialign_program: muscle-fast
