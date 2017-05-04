=======
Changes
=======

v0.7.0 (2017-05-04)
-------------------

* Add an “allele ratio” criterion to the germline filter to further reduce
  the number of false positives. The filter is activated by default and can
  be configured through the ``allele_ratio`` setting in the configuration
  file. :ref:`See the documentation for how it works <allele-ratio>.`
* Ignore the CDR3-encoding bases whenever comparing two V gene sequences.
* Avoid finding 5'-truncated V genes by extending found hits towards the
  5' end.
* By default, candidate sequences are no longer merged if they are nearly
  identical. That is, the ``differences`` setting within the two germline
  filter configuration sections is now set to zero by default.
  Previously, we believed the merging would remove some false
  positives, but it turns out we also miss true positives. It also seems
  that with the other changes in this version we also no longer get the
  particular false positives the setting was supposed to catch.
* Implement an experimental ``discoverj`` script for J gene discovery.
  It is curently not run automatically as part of ``igdiscover run``. See
  ``igdiscover discoverj --help`` for how to run it manually.
* Add a ``config`` subcommand, which can be used to change the
  configuration file from the command-line.
* Add a ``V_CDR3_start`` column to the ``assigned.tab``/``filtered.tab``
  tables. It describes where the CDR3 starts within the V sequence.
* Similarly, add a ``CDR3_start`` column to the ``new_V_germline.tab``
  file describing where the CDR3 starts within a discovered V sequence.
  It is computed by using the most common CDR3 start of the
  sequences within the cluster.
* Rename the ``compose`` subcommand to ``germlinefilter``.
* The ``init`` subcommand automatically fixes certain problems in the
  input database (duplicate sequences, empty records, duplicate sequence
  names). Previously, it would complain, but the user would have to fix
  the problems themselves.
* Move source code to GitHub
* Set up automatic code testing (continuous integration) via Travis
* Many documentation improvements

v0.6.0 (2016-12-07)
-------------------

* The FASTA files of the input V/D/J gene lists now need to be
  named ``V.fasta``, ``D.fasta`` and ``J.fasta``. The species name
  is no longer part of the file name. This should reduce confusion
  when working with species not supported by IgBLAST.
* The ``species:`` configuration setting in the configuration can
  (and should) now be left empty. Its only use was that it is passed
  to IgBLAST, but since IgDiscover provides IgBLAST with its own
  V/D/J sequences anyway, it does not seem to make a difference.
* A “cross-mapping” detection has been added, which should reduce
  the number of false positives.
  :ref:`See the documentation for an explanation <cross-mapping>`.
* Novel sequences identical to a database sequence no longer get the
  ``_S1234`` suffix.
* No longer trim trim the initial ``G`` run in sequences (due to RACE) by
  default. It is now a configuration setting.
* Add ``cdr3_location`` configuration setting: It allows to set whether to
  use a CDR3 in addition to the barcode for grouping sequences.
* Create a ``groups.tab.gz`` file by default (describing the de-barcoded
  groups)
* The pre-processing filter is now configurable. See the
  ``preprocessing_filter`` section in the configuration file.
* Many improvements to the documentation
* Extended and fixed unit tests. These are now run via a CI system.
* Statistics in JSON format are written to ``stats/stats.json``.
* IgBLAST 1.5.0 output can now be parsed. Parsing is also faster by 25%.
* More helpful warning message when no sequences were discovered in
  an iteration.
* Drop support for Python 3.3.

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
