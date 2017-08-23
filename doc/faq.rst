Questions and Answers
=====================


How many sequences are needed to discover germline V gene sequences?
--------------------------------------------------------------------

Library sizes of several hundred thousand sequences are required for V gene discovery, with even
higher numbers necessary for full database production. For example, IgM library sizes of 750,000
to 1,000,000 sequences for heavy chain databases and 1.5 to 2 million sequences for light chain
databases.


Can IgDiscover analyze IgG libraries?
-------------------------------------

IgDiscover has been developed to identify germline databases from libraries that contain
substantial fractions of unswitched antibody sequences. We recommend IgM libraries for heavy
chain V gene identification and IgKappa and IgLambda libraries for light chain identification.
IgDiscover can identify a proportion of gemline sequences in IgG libraries but the process is
much more efficient with IgM libraries, enabling the full set of germline sequences to be
discovered.


Can IgDiscover analyze a previously sequenced library?
------------------------------------------------------

Yes, IgDiscover accepts both unpaired FASTQ files and paired FASTA files but the program should
be made aware which is being used, see :ref:`input requirements <input-requirements>`.


Do the positions of the PCR primers make a difference to the output?
--------------------------------------------------------------------

Yes. For accurate V gene discovery, all primer sequences must be external to the V gene sequences.
For example, forward multiplex amplification primers should be present in the leader sequence or
5' UTR, and reverse amplification primers should be located in the constant region, preferably
close to the 5' border of the constant region. Primers that are present in framework 1 region or
J segments are not recommended for library production.


What are the advantages to 5'-RACE compared to multiplex PCR for IgDiscover analysis?
-------------------------------------------------------------------------------------

Both 5'-RACE and multiplex PCR have their own advantages.

5'-RACE will enable library production from species where the upstream V gene sequence is unknown.
The output of the ``upstream`` subcommand in IgDiscovery enables the identification of consensus
leader and 5'-UTR sequences for each of the identified germline V genes, that can subsequenctly
be used for primer design for either multiplex PCR or for monoclonal antibody amplification sets.

Multiplex PCR is recommended for species where the upstream sequences are well characterized.
Multiplex amplification products are shorter than 5'-RACE products and therefore will be easier
to pair and will have less length associated sequence errors.


What is meant by 'starting database'?
-------------------------------------

The starting database refers to the folder that contains the three FASTA files necessary for the
process of iterative V gene discovery to begin. IgDiscover uses the standalone IgBLAST program for
comparative assignment of sequences to the starting database. Because IgBlast requires three
files (for example ``V.fasta``, ``D.fasta``, ``J.fasta``), three FASTA files should be included
in the database folder for each analysis to proceed.

In the case of light chains (that do not contain D segments), a dummy D segment file should be
included as IgBLAST will not proceed if it does not see three files in the database folder. It is
sufficient to save the following sequence as a fasta file and rename it D.fasta, for example,
for it to function as the dummy ``D.fasta`` file for human light chain analysis::

    >D_ummy
    GGGGGGGGGG


How can I use the IMGT database as a starting database?
-------------------------------------------------------

Since we do not have permission to distribute IMGT database files with IgDiscover, you need to
download them directly from `IMGT <http://www.imgt.org/>`_.
See the :ref:`section about obtaining a V/D/J database <obtaining-database>`.


How do I change the parameters of the program?
----------------------------------------------

By editing :ref:`the configuration file <configuration>`.


Where do I find the individualized database produced by IgDiscover?
-------------------------------------------------------------------

The final germline database in FASTA format is in your :ref:`analysis
directory <analysis-directory>` in the subdirectory ``final/database/``. The ``V.fasta`` file
contains the new list of V genes. The ``D.fasta`` and ``J.fasta`` files are unchanged from the
starting database.

A phylogenetic tree of the V sequences can be found in ``final/dendrogram_V.pdf``.

For more details of how that database was created, you need to inspect the files created in the last
iteration of the discovery process, located in ``iteration-xx``, where ``xx`` is the number of
iterations configured in the ``igdiscover.yaml`` configuration file. For example, if three
iterations were used, look into ``iteration-03/``.

Most interesting in that folder are likely

- the linkage cluster analysis plots in ``iteration-03/clusterplots/``,
- the error histograms in ``iteration-03/errorhistograms.pdf``, which contain the windowed cluster
  analysis figures.
- Details about the individualized database in ``new_V_germline.tab`` in tab-separated-value format

The ``new_V_germline.fasta`` file is identical to the one in ``final/database/V.fasta``


What does the _S1234 at the end of same gene names mean?
--------------------------------------------------------

Please see the :ref:`Section on gene names <gene-names>`.