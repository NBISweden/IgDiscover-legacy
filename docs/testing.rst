.. _test:

=============
Test data set
=============

After installing IgDiscover, you should run it once on a small test data that we
provide, both to test your installation and to familiarize yourself with
running the program.

1. Download und unpack the test data set igdiscover-test-0.1.tar.gz_. To do this
   from the command-line, use these commands::

	wget https://export.uppmax.uu.se/b2014342/bc795179e13e84ef/igdiscover-test-0.1.tar.gz
	tar xf igdiscover-test-0.1.tar.gz

.. _igdiscover-test-0.1.tar.gz: https://export.uppmax.uu.se/b2014342/bc795179e13e84ef/igdiscover-test-0.1.tar.gz

   The test data set contains some paired-end reads from a rhesus monkey and the
   publicly available VH, DH, and JH gene databases as three FASTA files.

2. Initialize the IgDiscover pipeline directory::

	igdiscover init --db igdiscover-test-0.1/db/ --reads igdiscover-test-0.1/reads.1.fastq.gz discovertest

   The name ``discovertest`` is the name of the pipeline directory that will be
   created. Note that only the path to the *first* reads file needs to be
   given. The second file is found automatically.

   The command will have printed a message telling you that the pipeline
   directory has been initialized, that you should edit the configuration file,
   and how to actually run IgDiscover after that.

3. Edit the configuration file. The configuration file is in
   ``discovertest/igdiscover.yaml``. The configuration is in YAML format, which
   should be straightforward to understand. When editing the file, just follow
   the way it is already structured.

   Change the following settings in the configuration file:

   - Change the number of iterations from 4 to 2.
   - Change the ``barcode_length`` from 0 to 12.
   - Specify that a stranded protocol is used by changing the appropriate
     setting from *false* to *true*.
   - In the list of forward_primers, make sure that it contains this primer
     sequence: ``CGTGAGCTGAGTACGACTCACTATAGCTTCAC`` (The dash in the beginning
     that the other entries have is necessary. It marks each list item.)
   - In the same way, make sure that the ``reverse_primers`` contain
     ``GCAGGCCTTTTTGGCCNNNNNGGGGCATTCTCACAGGAGACGAGGGGGAAAAG``.

4. Run the analysis by changing into the pipeline directory and runing this command::

	cd discovertest && igdiscover run

5. Finally, inspect the results in the ``discovertest/final`` directory.
   See the :ref:`explanation of final result files <final-results>`.
