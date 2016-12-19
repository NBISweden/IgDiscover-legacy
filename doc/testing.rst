.. _test:

=============
Test data set
=============

After installing IgDiscover, you should run it once on a small test data that we
provide, both to test your installation and to familiarize yourself with
running the program.

1. Download und unpack `the test data set (version 0.3)`_. To do this
   from the command-line, use these commands::

	wget https://bitbucket.org/igdiscover/testdata/get/0.3.tar.gz
	tar xvf 0.3.tar.gz

.. _the test data set (version 0.3): https://bitbucket.org/igdiscover/testdata/get/0.3.tar.gz

   The test data set contains some paired-end reads from a rhesus monkey and a
   simplified version of the publicly available VH, DH, and JH gene lists
   as three FASTA files.

2. Initialize the IgDiscover pipeline directory::

       igdiscover init --db igdiscover-testdata-*/db/ --reads igdiscover-testdata-*/reads.1.fastq.gz discovertest

   The name ``discovertest`` is the name of the pipeline directory that will be
   created. Note that only the path to the *first* reads file needs to be
   given. The second file is found automatically.

   The command will have printed a message telling you that the pipeline
   directory has been initialized, that you should edit the configuration file,
   and how to actually run IgDiscover after that.

3. Instead of editing the configuration by hand, copy the prepared configuration
   file that comes with the test dataset into the pipeline directory::

       cp igdiscover-testdata-*/igdiscover.yaml discovertest/

   For your own runs, you will need to read through the configuration file
   and adjust it to your needs. You may want to have a look at it when the
   pipeline is running in the next step. The configuration is in YAML
   format. When editing the file, just follow the way it is already structured.

4. Run the analysis. To do so, change into the pipeline directory and run this
   command::

	cd discovertest && igdiscover run

   On this small dataset, running the pipeline should take not more than about 5 minutes.

5. Finally, inspect the results in the ``discovertest/final`` directory.
   For example, the final list of discovered V genes is in ``discovertest/final/database/V.fasta``.

   See the :ref:`explanation of final result files <final-results>`.
