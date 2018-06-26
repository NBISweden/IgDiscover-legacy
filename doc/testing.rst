.. _test:

=============
Test data set
=============

After installing IgDiscover, you should run it once on a small test data that we
provide, both to test your installation and to familiarize yourself with
running the program.

1. Download und unpack `the test data set (version 0.5)`_. To do this
   from the command-line, use these commands::

	wget https://bitbucket.org/igdiscover/testdata/downloads/igdiscover-testdata-0.5.tar.gz
	tar xvf igdiscover-testdata-0.5.tar.gz

.. _the test data set (version 0.5): https://bitbucket.org/igdiscover/testdata/downloads/igdiscover-testdata-0.5.tar.gz

   The test data set contains some paired-end reads from human IgM heavy chain
   dataset ERR1760498 and a database of IGHV, IGHD, IGHJ sequences based on
   Ensembl annotations. You should use a database of higher quality for your
   own experiments.

2. Initialize the IgDiscover pipeline directory::

       igdiscover init --db igdiscover-testdata/database/ --reads igdiscover-testdata/reads.1.fastq.gz discovertest

   The name ``discovertest`` is the name of the pipeline directory that will be
   created. Note that only the path to the *first* reads file needs to be
   given. The second file is found automatically. There may be a couple of
   messages “Skipping 'x' because it contains the same sequence as 'y'”, which
   you can ignore.

   The command will have printed a message telling you that the pipeline
   directory has been initialized, that you should edit the configuration file,
   and how to actually run IgDiscover after that.

3. The generated ``igdiscover.yaml`` configuration file does not actually need
   to be edited for the test dataset, but you may still want to have a read
   through it as you will need to do so for you own data. You may want to do
   this while the pipeline is running in the next step. The configuration is in
   YAML format. When editing the file, just follow the way it is already
   structured.

4. Run the analysis. To do so, change into the pipeline directory and run this
   command::

	cd discovertest && igdiscover run

   On this small dataset, running the pipeline should take not more than about 5 minutes.

5. Finally, inspect the results in the ``discovertest/iteration-01`` or
   ``discovertest/final`` directories. The discovered V genes and extra
   information are listed in
   ``discovertest/iteration-01/new_V_germline.tab``. Discovered J genes are
   in ``discovertest/iteration-01/new_J.tab``. There are also corresponding
   ``.fasta`` files with the sequences only.

   See the :ref:`explanation of final result files <final-results>`.


Other test data sets
--------------------

ENA project `PRJEB15295 <https://www.ebi.ac.uk/ena/data/view/PRJEB15295>`_ contains the data for
our Nature Communications paper from 2016, in particular
`ERR1760498 <https://www.ebi.ac.uk/ena/data/view/ERR1760498>`_, which is the data for the human “H1”
sample (multiplex PCR, IgM heavy chain).

Data used for testing TCR detection (human, RACE): `SRR2905677 <https://www.ncbi.nlm.nih.gov/sra/SRR2905677/>`_ and
`SRR2905710 <https://www.ncbi.nlm.nih.gov/sra/SRR2905710/>`_.
