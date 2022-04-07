=======
Changes
=======

v0.15 (2022-04-07)
------------------

* Change the algorithm used for describing how the discovered V gene differs from the
  germline gene (the ``database_changes`` column). This gives more sensible descriptions when
  the V gene is truncated at one end.
* Faster startup time (mostly noticable when using ``--version`` or ``--help``)
* Ensure candidates get a unique name even if the hashes (``_Sxxxx``) collide
* #108 Print a sensible error message when the GUI cannot be started.

v0.14 (2022-03-10)
------------------

* Fix a crash (``KeyError``) during "igdiscover augment" when region info
  for a database sequence could not be obtained.


v0.13 (2022-02-21)
------------------

* IgDiscover now uses AIRR-formatted files:
  See the `AIRR rearrangement schema <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`_
* IgBLAST is run with the appropriate parameters to produce AIRR-compliant files
* ``assigned.tab.gz`` and ``filtered.tab.gz`` contain this IgBLAST output plus extra columns
  that IgDiscover needs (the AIRR schema allows extra columns)
*``assigned.tab.gz`` and ``filtered.tab.gz`` are now called ``assigned.tsv.gz`` and
  ``filtered.tsv.gz`` (the ``.tsv`` extension is required by the AIRR specification)
* One downside is that, because there are more columns than before, the "assigned" and "filtered"
  files are larger than before.
* The upside is that these files can be used with other tools that accept AIRR-compliant files.
* Old "assigned" and "filtered" files can still be read by most IgDiscover commands. Output will
  always use new column names.
* The ``VDJ_nt`` column was removed to reduce file size somewhat. It is now recomputed when
  necessary from the appropriate offsets.
* Update to IgBLAST 1.17


v0.12 (2020-01-20)
------------------

* The ``discoverj`` command was renamed to ``discoverjd`` to reflect that it
  also supports D gene discovery.
* Previously, the ``why_filtered`` column would show a generic ``is_duplicate``
  reason for filters that compare candidates to each other. Now each filter
  criterion can be distinguised.
* The somewhat vague “too similar sequence” germline filter criterion
  incorrectly removed some candidates that have a mutation close to the 3’ end.
  This was replaced with a simpler filter that only ensures that there are no
  two candidates with the same sequence.
* Use IgBLAST 1.10
* Get rid of some unnecessary dependencies by no longer requiring the
  unmaintained ``sqt`` library. Installation with Conda is now faster and
  requires half the disk space.
* Add a *full_exact* column to ``candidates.tab``


v0.11 (2018-11-27)
------------------

* The IgBLAST cache is now disabled by default. We assume that, in most cases,
  datasets will not be re-run with the exact same parameters, and then it only
  fills up the disk. Delete your cache with ``rm -r ~/.cache/igdiscover`` to
  reclaim the space. To enable the cache, create a file
  ``~/.config/igdiscover.conf`` with the contents ``use_cache: true``.
* If you choose to enable the cache, results from the PEAR merging step will
  now also be cached. See also the :ref:`caching documentation <caching>`.
* Added detection of chimeras to the (pre-)germline filters. Any novel allele
  that can be explained as a chimera of two unmodified reference alleles is
  marked in the ``new_V_germline.tab`` file. This is a bit sensitive, so the
  candidate is currently not discarded.
* Two additional files ``annotated_V_germline.tab`` and
  ``annotated_V_pregermline.tab`` are created in each iteration during the
  germline filtering step. These are identical to the ``candidates.tab``
  file, except that they contain a ``why_filtered`` column that describes
  why a sequence was filtered. See the :ref:`documentation for this feature
  <annotated_v_tab>`.
* A more realistic test dataset (v0.5), now based on human instead of rhesus
  data, was prepared. The :ref:`testing instructions <test>` have been
  updated accordingly.
* J discovery has been tuned to give fewer truncated sequences.
* Statistics are written to ``stats/stats.json``.
* V SHM distribution plots are created automatically and written written to
  ``v-shm-distributions.pdf`` in each iteration folder.
* An ``igdiscover dbdiff`` subcommand was added that can compare two FASTA
  files.


v0.10 (2018-05-11)
------------------

* When computing a consensus sequence, allow some sequences to be truncated in
  the 3' end. Many of the discovered novel V alleles were truncated by one
  nucleotide in the 3' end because IgBLAST does not always extend the
  alignment to the end of the V sequence. If these slightly too short V
  sequences were in the majority, their consensus would lead to a truncated
  sequence as well. The new consensus algorithm allows for this effect at the
  3' end and can therefore more often than previously find the full sequence.
  Example::

     TACTGTGCGAGAGA (seq 1)
     TACTGTGCGAGAGA (seq 2)
     TACTGTGCGAGAG- (seq 3)
     TACTGTGCGAG--- (seq 4)
     TACTGTGCGAG--- (seq 5)

     TACTGTGCGAGAG  (previous consensus)
     TACTGTGCGAGAGA (new consensus)
* Add a column ``database_changes`` to the ``new_V_germline.tab`` file that
  describes how the novel sequence differs from the database sequence. Example:
  ``93C>T; 114A>G``
* Allow filtering by ``CDR3_shared_ratio`` and do so by default (needs
  documentation)
* Cache the edit distance when computing the distance matrix. Speeds up the
  ``discover`` command slightly.
* ``discover``: Use more than six CPU cores if available
* ``igblast``: Print progress every minute

v0.9 (2018-03-22)
-----------------

* Implemented allele ratio filtering for J gene discovery
* J genes are discovered as part of the pipeline (previously, one needed
  to run the ``discoverj`` script manually)
* In each iteration, dendrograms are now created not only for V genes, but
  also for D and J genes. The file names are ``dendrogram_D.pdf``,
  ``dendrogram_J.pdf``
* The V dendrograms are now in ``dendrogram_V.pdf`` (no longer
  ``V_dendrogram.pdf``). This puts all the dendrograms together when looking
  at the files in the iteration directory.
* The ``V_usage.tab`` and ``V_usage.pdf`` files are no longer created.
  Instead, ``expressed_V.tab`` and ``expressed_V.pdf`` are created. These
  contain similar information, but an allele-ratio filter is used to
  filter out artifacts.
* Similarly, ``expressed_D.tab`` and ``expressed_J.tab`` and their
  ``.pdf`` counterparts are created in each iteration.
* Removed ``parse`` subcommand (functionality is in the ``igblast`` subcommand)
* New CDR3 detection method (only heavy chain sequences): CDR3 start/end coordinates
  are pre-computed using the database V and J sequences. Increases detection rate
  to 99% (previously less than 90%).
* Remove the ability to check discovered genes for required motifs. This has never
  worked well.
* Add a column ``clonotypes`` to the ``candidates.tab`` that tries to count how many
  clonotypes are associated with a single candidate (using only exact occurrences).
  This is intended to replace the ``CDR3s_exact`` column.
* Add an ``exact_ratio`` to the germline filtering options. This checks the ratio
  between the exact V occurrence counts (``exact`` column) between alleles.
* Germline filtering option ``allele_ratio`` was renamed to ``clonotypes_ratio``
* Implement a cache for IgBLAST results. When the same dataset is re-analyzed,
  possibly with different parameters, the cached results are used instead of
  re-running IgBLAST, which saves a lot of time. If the V/D/J database or the
  IgBLAST version has changed, results are not re-used.

v0.8.0 (2017-06-20)
-------------------

* Add a ``barcodes_exact`` column to the candidates table. It gives the number
  of unique barcode sequences that were used by the sequences in the set of
  exact sequences. Also, add a configuration setting ``barcode_consensus``
  that can turn off consensus taking of barcode groups, which needs to be
  set to ``false`` for ``barcodes_exact`` to work.
* Add a ``Ds_exact`` column to candidates table.
* Add a ``D_coverage`` configuration option.
* The pre-processing filtering step no longer reads in the full table of
  IgBLAST assignments, but filters the table piece by piece. Memory usage
  for this step therefore does not depend anymore on the dataset size and
  should always be below 1 GB.
* The functionality of the ``parse`` subcommand has been integrated into
  the ``igblast`` subcommand. This means that ``igdiscover igblast`` now
  directly outputs a result table (``assigned.tab``). This makes it easier
  to use that subcommand directly instead of only via the workflow.
* The ``igblast`` subcommand now always runs ``makeblastdb`` by itself
  and deletes the BLAST database afterwards. This reduces clutter and
  ensures the database is always up to date.
* Remove the ``library_name`` configuration setting. Instead, the
  ``library_name`` is now always the same as the name of analysis
  directory.

v0.7.0 (2017-05-04)
-------------------

* Add an “allele ratio” criterion to the germline filter to further reduce
  the number of false positives. The filter is activated by default and can
  be configured through the ``allele_ratio`` setting in the configuration
  file. :ref:`See the documentation for how it works <allele-ratio>`.
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
