=======
Changes
=======

development version
-------------------
* Add column ``has_stop`` to ``candidates.tab``, which indicates whether the
  candidate sequence contains a stop codon.
* Add a configuration option that makes it possible to disable the 5' motif
  check by setting ``check_motifs: false`` (the ``looks_like_V`` column is
  ignored in this case).
* It should now be possible to install and run IgDiscover on OS X. This did
  not require changes in IgDiscover itself, but some dependencies needed to
  be made available as Conda packages (IgBLAST in particular).
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

v0.3
----
* Optionally, sequences are not renamed in the ``assigned.tab`` file, but
  retain their original name as in the FASTA or FASTQ file. Set ``rename:
  false`` in the configuration file to get this behavior.
* Started an “advanced” section in the manual.

v0.2
----

* IgDiscover can now also detect kappa and lambda light chain V genes (VK, VL)
