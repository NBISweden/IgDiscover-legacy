=======
Changes
=======

v0.5 (2016-09-01)
-----------------

* V sequences of the input database are now whitelisted by default.
  The meaning of the ``whitelist`` configuration option has changed:
  If set to ``false``, those sequences are no longer whitelisted.
  To whitelist additional sequences, create a ``whitelist.fasta``
  file as before.
* Sequences with stop codons are now filtered out by default.
* Use more stringent germline filtering parameters by default.

v0.4 (2016-08-24)
-----------------

* It is now possible to install and run IgDiscover on OS X. Appropriate Conda
  packages are available on bioconda.
* Add column ``has_stop`` to ``candidates.tab``, which indicates whether the
  candidate sequence contains a stop codon.
* Add a configuration option that makes it possible to disable the 5' motif
  check by setting ``check_motifs: false`` (the ``looks_like_V`` column is
  ignored in this case).
* Make it possible to whitelist known sequences: If a found gene candidate
  appears in that list, the sequence is included in the list of discovered
  sequences even when it would otherwise not pass filtering criteria. To enable
  this, just add a ``whitelist.fasta`` file to the project directory before
  starting the analysis.
* The criteria for germline filter and pre-germline filter are now configurable:
  See ``germline_filter`` and ``pre_germline_filter`` sections in the
  configuration file.
* Different runs of IgDiscover with the same parameters on the same input files
  will now give the same results. See the ``seed`` parameter in the configuration,
  also on how to get non-reproducible results as before.
* Both the germline and pre-germline filter are now applied in each iteration.
  Instead of the ``new_V_database.fasta`` file, two files named
  ``new_V_germline.fasta`` and ``new_V_pregermline.fasta`` are created.
* The ``compose`` subcommand now outputs a filtered version of the
  ``candidates.tab`` file in addition to a FASTA file. The table
  contains columns **closest_whitelist**, which is the name of the closest
  whitelist sequence, and **whitelist_diff**, which is the number of differences
  to that whitelist sequence.

v0.3 (2016-08-08)
-----------------

* Optionally, sequences are not renamed in the ``assigned.tab`` file, but
  retain their original name as in the FASTA or FASTQ file. Set ``rename:
  false`` in the configuration file to get this behavior.
* Started an “advanced” section in the manual.

v0.2
----

* IgDiscover can now also detect kappa and lambda light chain V genes (VK, VL)
