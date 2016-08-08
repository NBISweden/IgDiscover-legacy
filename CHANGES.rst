=======
Changes
=======

development version
-------------------
* Add column ``has_stop`` to ``candidates.tab``which indicates whether the
  candidate sequence contains a stop codon.
* Add a configuration option that makes it possible to disable the 5' motif
  check by setting ``check_motifs: false`` (the ``looks_like_V`` column is
  ignored in this case).

v0.3
----
* Optionally, sequences are not renamed in the ``assigned.tab`` file, but
  retain their original name as in the FASTA or FASTQ file. Set ``rename:
  false`` in the configuration file to get this behavior.
* Started an “advanced” section in the manual.

v0.2
----

* IgDiscover can now also detect kappa and lambda light chain V genes (VK, VL)
