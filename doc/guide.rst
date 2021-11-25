==========
User guide
==========


Overview
========

IgDiscover works on a single library at a time. It works within an
“analysis directory” for the library, which contains all intermediate
and result files.

To start an analysis, you need:

1. A FASTA or FASTQ file with single-end reads or two FASTQ files with
   paired-end reads (also, the files must be gzip-compressed)
2. A database of V/D/J genes (three FASTA files named ``V.fasta``, ``D.fasta``, ``J.fasta``)
3. A configuration file that describes the library

If you do not have a V/D/J database, yet, you may want to read the section about
:ref:`how to obtain V/D/J sequences <obtaining-database>`.

To run an analysis, proceed as follows.

.. note::
  If you are on macOS, it may be necessary to run ``export SHELL=/bin/bash`` before continuing.

1. Create and initialize the analysis directory.

   First, pick a name for your analysis. We will use ``myexperiment`` in the following.
   Then run ``igdiscover init``::

       igdiscover init myexperiment

   A dialog will appear and ask for the file with the *first* (forward) reads.
   Find your compressed FASTQ file that contains them and select it.
   Typical file names may be ``Library1_S1_L001_R1_001.fastq.gz`` or ``mylibrary.1.fastq.gz``.
   You do not need to choose the second read file!
   It is found automatically.

   Next, choose the directory with your database.
   The directory must contain the three files ``V.fasta``, ``D.fasta``, ``J.fasta``.
   These files contain the V, D, J gene sequences, respectively.
   Even if you have only light chains in your data, a ``D.fasta`` file needs to be provided;
   just use one with the heavy chain D gene sequences.

   If you do not want a graphical user interface, use the two command-line
   parameters ``--db`` and ``--reads1`` to provide this information instead::

       igdiscover init --db path/to/my/database/ --reads1 mylibrary.1.fastq.gz myexperiment

   Again, the second reads file will be found automatically.
   Use ``--single-reads`` instead of ``--reads1`` if you have single-end reads or a dataset with already merged reads.
   For ``--single-reads``, a FASTA file (not only FASTQ) is also allowed.
   In any case, an analysis directory named ``myexperiment`` will have been created.

2. Adjust the configuration file

   The previous step created a configuration file named ``myexperiment/igdiscover.yaml``, which
   you may :ref:`need to adjust <configuration>`.

3. Run the analysis

   Change into the newly created analysis directory and run the analysis::

       igdiscover run

   Depending on the size of your library, this will usually take a couple of hours. See the
   :ref:`running IgDiscover <running>` section for more fine-grained control over what to run and
   how to resume the process if something failed.


.. _obtaining-database:

Obtaining a V/D/J database
==========================

We use the term “database” to refer to three FASTA files that contain the sequences for the V, D
and J genes.
IMGT provides `sequences for download <http://www.imgt.org/vquest/refseqh.html>`_.
For discovering new VH genes, for example, you need to get the IGHV, IGHD and IGHJ files of your species.
As IgDiscover uses this only as a starting point, using a similar species will also work.

When using an IMGT database, it is very important to change the long IMGT sequence headers to
short headers as IgBLAST does not accept the long headers. You can use the program
``edit_imgt_file.pl``. If you installed IgDiscover from Conda, the script is already installed and
you can run it by typing the name. It is also
`available on the IgBlast FTP site <ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/>`_.

Run it for all three downloaded files, and then rename files appropriately to make sure that they
are named ``V.fasta``, ``D.fasta`` and ``J.fasta``.

You always need a file with D genes even if you analyze light chains.

In case you have used IgBLAST previously, note that there is *no need* to run the ``makeblastdb``
tool as IgDiscover will do that for you.


.. _input-requirements:

Input data requirements
=======================

Paired-end or single-end data
-----------------------------

IgDiscover can process input data of three different types:

* Paired-end reads in gzipped FASTQ format,
* Single-end reads in gzipped FASTQ format,
* Single-end reads in gzipped FASTA format.

IgDiscover was tested mainly on paired-end Illumina MiSeq reads (2x300bp), but it can also handle
454 and Ion Torrent data.

Depending on the input file type, use a variant of one of the following commands to initialize
the analysis directory::

    igdiscover init --single-reads=file.fasta.gz  --database=my-database-dir/ myexperiment
    igdiscover init --reads1=file.1.fasta.gz  --database=my-database-dir/ myexperiment
    igdiscover init --reads1=file.1.fastq.gz  --database=my-database-dir/ myexperiment


Read layout
-----------

Paired-end reads are first merged and then processed in the same way as single-end reads. Reads
that could not be merged are discarded. Single-end reads and merged paired-end reads are expected
to follow this structure (from 5' to 3'):

* The forward primer sequence. This is optional.
* A UMI (random barcode). This is optional. Set the configuration option ``barcode_length_5p`` to 0
  if you don’t have random barcodes or if you don’t want the program to use them.
* Optionally, a run of G nucleotides. This is an artifact of the RACE protocol (Rapid
  amplification of cDNA ends). If you have this, set ``race_g`` to ``true`` in the configuration file.
* 5' UTR
* Leader
* Re-arranged V, D and J gene sequences for heavy chains; only V and J for light chains
* An optional UMI (random barcode). Set the configuration option ``barcode_length_3p`` to the
  length of this UMI. You can currently not have both a 5' and a 3' UMI.
* The reverse primer. This is optional.

We use IgBLAST to detect the location of the V, D, J genes through the
``igdiscover igblast`` subcommand. The G nucleotides
after the barcode are split off if the configuration specifies
``race_g: true``. The leader sequence is detected by looking for a start
codon near 60 bp upstream of the start of the V gene match.



.. _configuration:

Configuration
=============

The ``igdiscover init`` command creates a configuration file
``igdiscover.yaml`` in the analysis directory. To configure
your analysis, change that file with a text editor before
running the analysis with ``igdiscover run``.


The syntax should be mostly self-explanatory.
The file is in YAML format, but you will not need to learn that.
Just follow the examples given in the file.
A few rules that may be good to know are the following ones:

1. Lines starting with the ``#`` symbol are comments (they are ignored)
2. A configuration option that is meant to be switched on or off will say something like ``stranded: false`` if it is off.
   Change this to ``stranded: true`` to switch the option on (and vice versa).
3. The primer sequences are given as a list, and must be written in a certain way - one sequence per line, and a ``-`` (dash) in front, like so::

       forward_primers:
       - ACGTACGTACGT
       - AACCGGTTAACC

   Even if you have only one primer sequence, you still need to use this syntax.

To find out what the configuration options achieve, see the explanations in the configuration file itself.

The main parameters that may require adjusting are the following.

The ``iterations`` option sets the number of rounds of V gene discovery
that will be performed. By default, one iteration is run. In each
iteration, all the sequences will be mapped with IgBLAST, which is the
most time-consuming part of running the pipeline. Thus, when you go from 1 to 2
iterations, you almost double the runtime requirements.

In previous IgDiscover versions, more iterations than one were necessary,
but we have improved sensitivity since then, so you should not need to increase
this.

Especially for nearly complete starting databases, for example when
analyzing a human IgM library with the current IMGT heavy chain database,
a single iteration is totally sufficient to produce an individualized database.

If you start with a very small V database (for example with only a single
V gene sequence), you may get better results when you increase this to 2.

If you do not want to discover any new genes, then use ``iterations: 0``.
This may be useful to only produce an expression profile, for example.

The ``ignore_j`` option should be set to ``true`` when producing a V gene
database for a species where J sequences are unknown::

    ignore_j: true

Setting the parameters ``stranded``, ``forward_primers`` and ``reverse_primers``
to the correct values can be used to remove 5' and 3' primers from the sequences.
Doing this is not strictly necessary for IgDiscover. It is simplest
if you do not specify any primer sequences.


Pregermline and germline filter criteria
----------------------------------------

IgDiscover V gene discovery works in two stages: The program first generates
a list of *candidate* V gene sequences. This list includes many false
positives. In the subsequent *germline filtering* step, the list is
therefore trimmed rigorously in order to produce the final list of germline
sequences. See also the `germline filtering <germline-filters>`:ref: section.

The stringency requirements for the germline filter can be set in the
configuration file in the `germline_filter` and `pregermline_filter`
sections.

The `pregermline_filter` section is used in all but the last iteration.
That is, it is ignored if you use the default of running only a single
iteration.

The idea behind the pregermline filter is to initially use less stringent
riteria, which allows to grow the starting database more quickly, but at
the risk of adding some false positives. The last iteration, in which the
more stringent `germline_filter` settings are used, will then remove those
remaining false positives.

Here is how it looks in the configuration file::

   pre_germline_filter:
     unique_cdr3s: 2      # Minimum number of unique CDR3s (within exact matches)
     unique_js: 2         # Minimum number of unique J genes (within exact matches)
     check_motifs: false  # Check whether 5' end starts with known motif
     whitelist: true      # Add database sequences to the whitelist
     cluster_size: 0      # Minimum number of sequences assigned to cluster
     allow_stop: true     # Whether to allow non-productive sequences containing stop codons
     cross_mapping_ratio: 0.02  # Threshold for removal of cross-mapping artifacts (set to 0 to disable)
     allele_ratio: 0.1    # Required minimum ratio between alleles of a single gene

   # Filtering criteria applied to candidate sequences in the last iteration.
   # These should be more strict than the pre_germline_filter criteria.
   #
   germline_filter:
     unique_cdr3s: 5      # Minimum number of unique CDR3s (within exact matches)
     unique_js: 3         # Minimum number of unique J genes (within exact matches)
     check_motifs: false  # Check whether 5' end starts with known motif
     whitelist: true      # Add database sequences to the whitelist
     cluster_size: 100    # Minimum number of sequences assigned to cluster
     allow_stop: false    # Whether to allow non-productive sequences containing stop codons
     cross_mapping_ratio: 0.02  # Threshold for removal of cross-mapping artifacts (set to 0 to disable)
     allele_ratio: 0.1    # Required minimum ratio between alleles of a single gene

Factors that affect germline discovery include library source (IgM vs IgK, IgL or IgG)
library size, sequence error rate and individual genomic factors (for example the
number of J segments present in an individual).

In general, setting a higher cutoff of ``unique_cdr3s`` and ``unique_js`` will minimize the number
of false positives in the output. Example::

   unique_cdr3s: 10      # Minimum number of unique CDR3s (within exact matches)
   unique_js: 4          # Minimum number of unique J genes (within exact matches)

Read also about the :ref:`cross mapping <cross-mapping>`, for which germline filtering corrects, and
about the :ref:`germline filters <germline-filters>`.

.. versionchanged::
   The ``differences`` configuration setting was removed.


.. _jdiscovery:

IgDiscover will also try to discover which J genes are used in the input sample. J discovery
is configured in the ``j_discovery`` section in the configuration file. It looks like this::

    j_discovery:
      allele_ratio: 0.2         # Required minimum ratio between alleles of a single gene
      cross_mapping_ratio: 0.1  # Threshold for removal of cross-mapping artifacts.
      propagate: true           # Use J genes discovered in iteration 1 in subsequent ones



.. _running:

Running IgDiscover
==================

Resuming failed runs
--------------------

The command ``igdiscover run``, which is used to start the pipeline, can also be used to resume
execution if there was an interruption (a transient failure). Reasons for interruptions might be:

* Ctrl+C was pressed on the keyboard
* A full harddisk
* If running on a cluster, the program may have been terminated because it exceeded its allocated
  time
* Too little RAM
* Power loss

To resume execution after you have fixed the problem, go to the analysis directory and run
``igdiscover run`` again. It will skip the steps that have already finished successfully.
This capability comes from the workflow management system
`snakemake <https://snakemake.bitbucket.io/>`_, on which ``igdiscover run`` is based.
Snakemake will determine automatically which steps need to be re-run in order to get to a full
result and then run only those.

Alterations to the configuration file after an interruption are possible, but affect only
steps that have not already finished successfully. For example, assume you interrupted a
run with Ctrl+C after it is already past the step in which barcodes are removed. Then,
even if you change the barcode length in the configuration, the barcode removal step will
not be re-run when you resume the pipeline and the previous barcode length is in effect.
See also the next section.

Changing parameters and re-running parts of the pipeline
--------------------------------------------------------

When you experiment with parameters in the ``igdiscover.yaml`` file, such as
germline filtering criteria, you do not need to re-run the entire pipeline from
the beginning, but can re-use the results that already exist. This can save a lot
of processing time, in particular when you avoid re-running IgBLAST in this way.

As described in the previous section, ``igdiscover run`` automatically figures out
which files need to be re-created if a run was interrupted. Unfortunately, this
mechanism is currently not smart enough to also look for changes in the
``igdiscover.yaml`` file. Thus, if the full pipeline has finished successfully,
then re-running ``igdiscover run`` will just print the message ``Nothing to be done.``
even after you have changed the configuration file.

You will therefore need to know yourself which file you want to regenerate.
Then follow the following steps. Note that these will remove parts of the existing
results, and if you need to keep them, make a copy of your analysis directory first.

1. Change the configuration setting.
2. Delete the file that needs to be re-generated. Assume it is ``filename``
3. Run ``igdiscover run filename`` to re-create the file. Only that file
   will be created, not the ones that usually would be created afterwards.
4. Optionally, run ``igdiscover run`` (without a file name this time) to
   update the remaining files (those that depend on the file that was just
   updated).

For example, assume you want to modify some germline filtering setting and then re-run
the pipeline. Change the setting in your ``igdiscover.yaml``, then run these
commands::

    rm iteration-01/new_V_germline.tab
    igdiscover run iteration-01/new_V_germline.tab

The above will have regenerated the ``iteration-01/new_V_germline.tab`` file
and also the ``iteration-01/new_V_germline.fasta`` file since they are
generated by the same script. If you want to update any other files, then also
run ::

    igdiscover run


.. _analysis-directory:

The analysis directory
======================

IgDiscover writes all intermediate files, the final V gene database, statistics and plots into
the analysis directory that was created with ``igdiscover init``.
Inside that directory, there is a ``final/`` subdirectory that contains the analysis results.

These are the files and subdirectories that can be found in the analysis directory.
Subdirectories are described in detail below.

igdiscover.yaml
    The configuration file.
    Make sure to adjust this to your needs as described above.

reads.1.fastq.gz, reads.2.fastq.gz
    Symbolic links to the raw paired-end reads.

database/
    The input V/D/J database (as three FASTA files).
    The files are a copy of the ones you selected when running ``igdiscover init``.

reads/
    Processed reads (merged, de-duplicated etc.)

iteration-xx/
    Iteration-specific analysis directory, where “xx” is a number starting from 01.
    Each iteration is run in one of these directories.
    The first iteration (in ``iteration-01``) uses the original input database, which is also found in the ``database/`` directory.
    The database is updated and then used as input for the next iteration.

final/
    After the last iteration, IgBLAST is run again on the input sequences, but using the final database (the one created in the very last iteration).
    This directory contains all the results, such as plots of the repertoire profiles.
    If you set the number of iterations to 0 in the configuration file, this directory is the only one that is created.


.. _final-results:

Final results
-------------

Final results are found in the ``final/`` subdirectory of the analysis directory.

final/database/(V,D,J).fasta
    These three files represent the final, individualized V/D/J database found by IgDiscover.
    The D file is a copy of the original starting database; it is not updated by IgDiscover.

final/dendrogram_(V,D,J).pdf
    These three PDF files contain dendrograms of the V, D and J sequences in the individualized
    database.

final/assigned.tsv.gz
    V/D/J gene assignments and other information for each sequence. This is an AIRR-compliant
    TSV file, obtained from running IgBLAST and adding some IgDiscovere-specific columns.
    See the `AIRR rearrangement specification <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`_.
    See also below for a description of the columns.

final/filtered.tsv.gz
    Filtered V/D/J gene assignments. This is the same as the ``assigned.tsv.gz``, but with
    low-quality assignments filtered out.
    Run ``igdiscover filter --help`` to see the filtering criteria.

final/expressed_(V,D,J).tab, final/expressed_(V,D,J).pdf
    The V, D and J gene expression counts. Some assignments are filtered out to reduce artifacts. In particular,
    an allele-ratio filter of 10% is applied. For D genes, only those with an E-value of at most
    1E-4 and a coverage of at least 70% are counted. See also the help for the ``igdiscover count``
    subcommand, which is used to create these files.

    The ``.tab`` file contains the counts as a table, while the PDF file contains a plot of the same values.

    These tables also exist in the iteration-specific directories (``iteration-xx``). For those,
    note that the numbers do not include the genes that were discovered in that iteration. For
    example, ``iteration-01/expressed_V.tab`` shows only expression counts of the V genes in the
    starting database.

final/errorhistograms.pdf
    A PDF with one page per V gene/allele.
    Each page shows a histogram of the percentage differences for that gene.

final/clusterplots/
    This is a directory that contains one PNG file for each discovered gene/allele.
    Each image shows a clusterplot of all the sequences assigned to that gene.
    Note that the shown clusterplots are by default restricted to showing only at most 300 sequences,
    while the actual clustering used by IgDiscover uses 1000 sequences.

If you are interested in the results of each iteration, you can inspect the iteration-xx/ directories.
They are structured in the same way as the final/ subdirectory, except that the results are based on the intermediate databases of that iteration.
They also contain the following additional files.

iteration-xx/candidates.tab
    A table with candidate novel V alleles (or genes).
    This is a list of sequences found through the *windowing strategy* or *linkage cluster analysis*,
    as discussed in our paper. See :ref:`the full description of candidates.tab <candidates_tab>`.

iteration-xx/read_names_map.tab
    For each candidate novel V allele listed in ``candidates.tab``, this file contains one row that
    lists which sequences went into generating this candidate. Only the exact matches are listed,
    that is, the number of listed sequence names should be equal to the value in the *exact*
    column. Each line in this file contains tab-separated values. The first is name of the
    candidate, the others are the names of the sequences. Some of these sequences may be consensus
    sequences if barcode grouping was enabled, so in that case, this will not be a read name.

iteration-xx/new_V_germline.fasta, iteration-xx/new_V_pregermline.fasta
    The discovered list of V genes for this iteration.
    The file is created from the ``candidates.tab`` file by applying either the germline or pre-germline filter.
    The file resulting from application of the germline filter is used in the last iteration only.
    The file resulting from application of the pre-germline filter is used in earlier iterations.

iteration-xx/annotated_V_germline.tab, iteration-xx/annotated_V_pregermline.tab
    A version of the ``candidates.tab`` file that is annotated with extra columns that describe why
    a candidate was filtered out. See :ref:`the description of this file <annotated_v_tab>`.

iteration-xx/new_J.tab, iteration-xx/new_J.fasta
    The discovered list of J genes for this iteration.


Other files
-----------

For completeness, here is a description of the files in the ``reads/`` and ``stats/`` directories.
They are created during pre-processing and are not iteration specific.

reads/1-limited.1.fastq.gz, reads/1-limited.1.fastq.gz
    Input reads file limited to the first N entries. This is just a symbolic
    link to the input file if the ``limit`` configuration option is not set.

reads/2-merged.fastq.gz
    Reads merged with PEAR or FLASH

reads/3-forward-primer-trimmed.fastq.gz
    Merged reads with 5' primer sequences removed. (This file is automatically removed when
    it is not needed anymore.)

reads/4-trimmed.fastq.gz
    Merged reads with 5' and 3' primer sequences removed.

reads/5-filtered.fasta
    Merged, primer-trimmed sequences converted to FASTA, and too short sequences removed.
    (This file is automatically removed when it is not needed anymore.)

reads/sequences.fasta.gz
    Fully pre-processed sequences. That is, filtered sequences without duplicates (using VSEARCH)

stats/reads.txt
    Statistics of pre-processed sequences.

stats/readlengths.txt, stats/readlengths.pdf
    Histogram of the lengths of pre-processed sequences (created from ``reads/sequences.fasta``)


Format of output files
======================


assigned.tsv.gz
---------------

This file is a gzip-compressed table with tab-separated values. It follows the
`AIRR Rearrangement Schema <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`_.
It is created by the ``igdiscover igblast`` subcommand followed by ``igdiscover augment``.
Columns ``sequence_id`` to ``np2_length`` have a meaning as per that schema and most of them are
copied unmodified from IgBLAST. Subsequent columns are specific to IgDiscover and are ignored
by other tools accepting AIRR-formatted tables.

In brief, the first row is a header with column names, and each subsequent row describes the IgBLAST
results for a single pre-processed input sequence.

Note: This file is typically quite large.
LibreOffice can open the file directly (even though it is compressed), but make sure you have enough RAM.

Columns specific to IgDiscover (added by ``igdiscover augment``):

count
    How many copies of input sequence this query sequence represents. Copied from the ``;size=3;``
    entry in the FASTA header field that is added by ``VSEARCH -derep_fulllength``.

V_covered, D_covered, J_covered
    percentage of bases of the reference gene that is covered by the bases of the query sequence

FR1_SHM, CDR1_SHM, FR2_SHM, CDR2_SHM, FR3_SHM, V_SHM, J_SHM
    rate of somatic hypermutation (actually, an error rate)

V_errors, J_errors
    Absolute number of errors (differences) in the V and J gene match

UTR
    Sequence of the 5' UTR (the part before the V gene match up to, but not including, the start codon)
    (Note: Currently not included.)

leader
    Leader sequence (the part between UTR and the V gene match)
    (Note: Currently not included.)

V_CDR3_start
    Start coordinate of CDR3 within ``V_nt``. Set to zero if no CDR3 was detected.
    Comparisons involving the V gene ignore those V bases that are part of the CDR3.

name, barcode, race_G, genomic_sequence
    see the following explanation

The UTR, leader, barcode, race_G and genomic_sequence columns are filled in the following way.

1. Split the 5' end barcode from the sequence (if barcode length is zero, this will be empty), put it in the **barcode** column.
2. Remove the initial run of G bases from the remaining sequence, put that in the **race_G** column.
3. The remainder is put into the **genomic_sequence** column.
4. If there is a V gene match, take the sequence *before* it and split it up in the following way. Search for the start codon and write the part before it into the **UTR** column. Write the part starting with the start column into the **leader** column.


filtered.tsv.gz
---------------

This table is the same as the ``assigned.tsv.gz`` table, except that rows containing low-quality matches have been filtered out.
Rows fulfilling any of the following criteria are filtered:

- The J gene was not assigned
- A stop was codon found
- The V gene coverage is less than 90%
- The J gene coverage is less than 60%
- The V gene E-value is greater than 10\ :sup:`-3`


.. _candidates_tab:

candidates.tab
--------------

This table contains the candidates for novel V genes found by the ``discover`` subcommand.
As the other files, it is a text file in tab-separated values format, with the first row containing the column headings.
It can be opened directly in LibreOffice, for example.

Candidates are found by inspecting all the sequences assigned to a database gene, and clustering them in multiple ways.
The candidate sequences are found by computing a consensus from each found cluster.

Each row describes a single candidate, but possibly multiple clusters.
If there are multiple clusters from a single gene that lead to the same consensus sequence, then they get only one row.
The *cluster* column lists the source clusters for the given sequence.
Duplicate sequences can still occur when two different genes lead to identical consensus sequences.
(These duplicated sequences are merged by the germline filters.)

Below, we use the term *cluster set* to refer to all the sequences that are in any of the listed clusters.

Some clusters lead to ambiguous consensus sequences (those that include ``N`` bases).
These have already been filtered out.


name
    The name of the candidate gene. See :ref:`novel gene names <gene-names>`.

source
    The original database gene to which the sequences from this row were originally assigned.
    All candidates coming from the same source gene are grouped together.

chain
    Chain type: *VH* for heavy, *VK* for light chain lambda, *VL* for light chain kappa

cluster
    From which type of cluster or clusters the consensus was computed.
    If there are multiple clusters that give rise to the same consensus sequence, they are all listed here, separated by semicolon.
    A cluster name such as ``2-4`` is for a percentage difference window:
    Such a cluster consists of all sequences assigned to the source gene that have a percentage difference to it between 2 and 4 percent.

    A cluster name such as ``cl3`` describes a cluster generated through linkage cluster analysis.
    The clusters are simply named ``cl1``, ``cl2``, ``cl3`` etc.
    If any cluster number seems to be missing (such as when cl1 and cl3 occur, but not cl2), then this means that the cluster led to an ambiguous consensus sequence that has been filtered out.
    Since the ``cl`` clusters are created from a random subsample of the data (in order to keep computation time down),
    they are never larger than the size of the subsample (currently 1000).

    The name ``db`` represents a cluster that is identical to the database sequence.
    If no actual cluster corresponding to the database sequence is found, but the database sequence is expressed, a ``db`` cluster is inserted artificially in order to make sure that the sequence is not lost.
    The cluster name ``all`` represents the set of all sequences assigned to the source gene.
    This means that an unambiguous consensus could be computed from all the sequences.
    Typically, this happens during later iterations when there are no more novel sequences among the sequences assigned to the database gene.

cluster_size
    The number of sequences from which the consensus was computed.
    Equivalently, the size of the cluster set (all clusters described in this row).
    Sequences that are in multiple clusters at the same time are counted only once.

Js
    The number of unique J genes associated with the sequences in the cluster set.

    Consensus sequences are computed only from V gene sequences, but each V gene sequence is part of a full V/D/J sequence.
    We therefore know for each V sequence which J gene it was found with.
    This number says how many different J genes were found for all sequences that the consensus in this row was computed from.

CDR3s
    The number of unique CDR3 sequences associated with the sequences in the cluster set.
    See also the description for the *Js* column.
    This number says how many different CDR3 sequences were found for all sequences that the consensus in this row was computed from.

exact
    The number of exact occurrences of the consensus sequence among all sequences assigned to the
    source gene, ignoring the 3' junction region.

    While the consensus sequence is computed only from a subset of sequences assigned to a source
    gene, *all* sequences assigned to the source gene are searched for exact occurrences
    of that consensus sequence.

    When comparing sequences, they are first truncated at the 3' end by removing those (typically
    8) bases that correspond to the CDR3 region.

full_exact
    The number of full, exact occurrences of the consensus among all sequences assigned to the
    source gene. This is the same as the *exact* column, but without removing the 3' junction
    region.

barcodes_exact
    How many unique barcode sequences were used by the sequences in the set of exact sequences
    (described above).

Ds_exact
    How many unique D genes were used by the sequences in the set of exact sequences (described
    above). Only those D gene assignments are included in this count for which the number of errors
    is zero, the E-value is at most a given threshold, and for which the number of covered bases
    is at least a given percentage.

Js_exact
    How many unique J genes were used by the sequences in the set of exact sequences (described above).

CDR3s_exact
    How many unique CDR3 sequences were used by the sequences in the set of exact sequences (described above).

clonotypes
    The estimated number of clonotypes within the set of exact sequences (which is described above).
    The value is computed by clustering the unique CDR3 sequences associated with all exact
    occurrences, allowing up to six differences (mismatches, insertions, deletions) and then
    counting the number of resulting clusters.

database_diff
    The number of differences between the consensus sequence and the sequence of the source gene.
    (Given as edit distance, that is insertion, deletion, mismatch count as one difference each.)

has_stop
    Indicates whether the consensus sequence contains a stop codon.

looks_like_V
    Whether the consensus sequence “looks like” a true V gene (1 if yes, 0 if no).
    Currently, this checks whether the 5' end of the sequence matches a known V gene motif.

CDR3_start
    Where the CDR3 starts within the discovered V gene sequence. This uses the most common
    CDR3 start location among the sequences from which this consensus is derived.

consensus
    The consensus sequence itself.

The ``igdiscover discover`` command can also be run by hand with other parameters, in which case additional columns may appear.

N_bases
    Number of ``N`` bases in the consensus


.. _annotated_v_tab:

annotated_V_*.tab
-----------------

The two files ``annotated_V_germline.tab`` and ``annotated_V_pregermline.tab`` are copies of the
``candidates.tab`` file with two extra columns that show *why* a candidate was filtered in the
germline and pre-germline filtering steps. The two columns are:

  * ``is_filtered`` – A number describing how many filtering criteria exclude this candidate.
  * ``why_filtered`` – A semicolon-separated list of filtering reasons.

The following values can occur in the ``why_filtered`` column:

too_low_dbdiff
    The number of differences between this candidate and the database is lower than the required number.

too_many_N_bases
    The candidate contains too many ``N`` wildcard characters.

too_low_CDR3s_exact
    The ``CDR3s_exact`` value for this candidate is lower than the configured threshold.

too_high_CDR3_shared_ratio
    The ``CDR3_shared_ratio`` is higher than the configured threshold.

too_low_Js_exact
    The ``Js_exact`` value is lower than the configured threshold.

has_stop
    The filter configuration disallows stop codons, but this candidate has one and is not whitelisted.

too_low_cluster_size
    The ``cluster_size`` of this candidate is lower than the configured threshold, and the candidate is not whitelisted.

xmap_ratio
    The :ref:`cross-mapping ratio <cross-mapping>` between this candidate and another one is too
    high. This is written as ``xmap_ratio=0.01,other=VH1-1``, where 0.01 is the cross-mapping ratio
    and “VH1-1” is the sequence to which the comparison was made.

clonotype_ratio
    The :ref:`clonotype ratio <allele-ratio>` between two alleles of the same gene is too low.
    This is written as ``clonotype_ratio=0.03,other=VH1-1``, where 0.03 is the clonotype ratio
    and “VH1-1” is the name of the sequence to which the comparison was made.

ex_occ_ratio
    The :ref:`exact occurrence ratio <allele-ratio>` between two alleles of the same gene is too
    low. This is written as ``ex_occ_ratio=0.03,other=VH1-1``, where 0.03 is the exact occurrence
    ratio and “VH1-1” is the name of the sequence to which the comparison was made.

Ds_exact_ratio
    The ratio of Ds_exact_ratio between two alleles of the same gene is too low. This is written as
    ``Ds_exact_ratio=0.03,other=VH1-1``, where 0.03 is the ratio and VH1-1 is the name of the
    sequence to which the comparison was made.

identical_to
    The candidate is identical to another one (a duplicate) or a truncated version of another one.
    This is written as ``identical_to=VH1-1,truncated``, where “VH1-1” is the name of the other
    sequence. If the ``truncated`` part is missing, then the sequences were exactly identical.


.. _stats-json:

stats/stats.json
----------------

The ``stats/stats.json`` is a JSON file that contains various statistics about a discovery run.
In particular, below the ``"iterations"`` and then the ``"database"`` key, you will see information
about the number of found V alleles found or lost in each iteration.

Note that both the germline filter and pregermline filter are run in each iteration, resulting
in two databases, and therefore the ``stats.json`` file contains information about both. This
allows one to compare the two filters. That is, the germline-filtering numbers
tell you what the size of the database would be *if this was the last iteration*, but if there
are more iterations, then the numbers for the pregermline-filtered database are the relevant ones.

Here is a shortened example::

    "iterations": [
      {
        "database": {
          "iteration": 0,
          "size": 44
        }
      },
      {
        "assignment_filtering": {
          (omitted)
        },
        "database": {
          "iteration": 1,
          "size": 12,
          "gained": 7,
          "lost": 39,
          "size_pre": 20,
          "gained_pre": 15,
          "lost_pre": 39
        }
      }
    ]

iteration
    The iteration number. The first entry is iteration zero and is not an actual iteration, but
    gives information about the size of the starting database.

size
    The number of V alleles in the germline-filtered database that was discovered in this iteration.
    For iteration zero, this gives the number of V alleles in the starting database.

size_pre
    Same as ``size``, but for the pregermline-filtered version of the database.

gained
    Number of novel V alleles in the germline-filtered database that was discovered in this
    iteration, compared to the pregermline-filtered database of the previous iteration.

gained_pre
    Same as ``gained``, but the comparison is made between the current and previous
    pregermline-filtered databases.

lost
    The number of alleles that existed in the pregermline-filtered database of the previous
    iteration, but are not present in the current germline-filtered database.

lost_pre
    Same as ``lost``, but the comparison is made between the current and previous
    pregermline-filtered databases.

Let us look at the above example.

* ``size``=12: If this were the last iteration, IgDiscover would give a
  final database with 12 V alleles.
* ``gained``=7: 7 of those 12 alleles are novel compared to the database
  of the previous iteration (the starting database in this case).
* ``lost``=39: Of the 44 alleles that were in the database of the previous iteration, 39 could
  not be found in this iteration (using the germline filter), that is, 5 alleles are common.
* ``size_pre``=20: If this is not the final iteration, IgDiscover starts the next iteration
  with an intermediate input database containing 20 alleles.
* ``gained_pre``=15: Of those 20 alleles, 15 are new compared to the previous iteration.

Other notes:

* The most important values are ``size`` and ``gained``.
* All four “gained” and “lost” values show comparisons to the pregermline-filtered database of the
  previous iteration.
* In iteration 1, when there is no previous iteration, the comparison is made to the
  starting database.


.. _gene-names:

Names for discovered genes
--------------------------

Each gene discovered by IgDiscover gets a unique name such as “VH4.11_S1234”.
The “VH4.11” is the name of the database gene to which the novel
V gene was initially assigned. The number *1234* (hash) is derived from the nucleotide
sequence of the novel gene. That is, if you discover the same sequence in two
different runs of the IgDiscover, or just in different iterations, the number will
be the same. This may help when manually inspecting results.

Be aware that you still need to check the sequence itself since even different
sequences can sometimes lead to the same number (a “hash collision”).

The ``_S1234`` suffixes do not accumulate.
Before IgDiscover adds the suffix in an iteration, it removes the suffix if it already exists.


Subcommands
===========

The ``igdiscover`` program has multiple subcommands.
You should already be familiar with the two commands ``init`` and ``run``.
Each subcommand comes with its own help page that shows how to use that subcommand.
Run the command with the ``--help`` option to see the help. For example, ::

    igdiscover run --help

shows the help for the ``run`` subcommand.

The following additional subcommands may be useful for further analysis.

commonv
    Find common V genes between two different antibody libraries

upstream
    Cluster upstream sequences (UTR and leader) for each gene

dendrogram
    Draw a dendrogram of sequences in a FASTA file.

rename
    Rename sequences in a target FASTA file using a template FASTA file

union
    Compute union of sequences in multiple FASTA files

:ref:`clonotypes <clonotypes>`
    List the clonotypes (unique V, J, CDR3 combinations) present in a sample


The following subcommands are used internally, and listed here for completeness.

filter
    Filter a table with IgBLAST results

count
    Count and plot V, D, J gene usage

group
    Group sequences by barcode and V/J assignment and print each group’s consensus (unused in IgDiscover)

germlinefilter
    Create new V gene database from V gene candidates using the germline and pre-germline filter
    criteria.

discover
    Discover candidate new V genes within a single antibody library

clusterplot
    For each V gene, plot a clustermap of the sequences assigned to it

errorplot
    Plot histograms of differences to reference V gene


.. _clonotypes:

``igdiscover clonotypes``
-------------------------

The ``igdiscover clonotypes`` command lists the clonotypes present in a sample.
The only required parameter is the name of a file with assigned sequences.
Normally, this will be a ``filtered.tsv.gz`` file.

Two sequences are considered to be of the same clonotype if

- their V and J assignments are the same
- the length of their CDR3 is identical
- their CDR3 sequences are similar (see below for what this means)

That is, clonotypes are found by clustering the input sequences by V gene,
J gene and CDR3 similarity (using single-linkage clustering).

For each cluster, a representative row (assignment) is chosen and
considered to be the clonotype. The output is a table with one row
per clonotype. It is written to standard output.

By default, the output table is sorted by V/D/J gene names.
Use ``--sort`` to sort by group size (largest first).

Similarity
~~~~~~~~~~

To determine whether two CDR3s are similar, the Hamming distance
between the CDR3 nucleotide sequences (``CDR3_nt`` column) must be
at most 1. To allow more differences, use ``--mismatches``. To
compare amino acid sequences (``CDR3_aa``) instead, use ``--aa``.

The members file
~~~~~~~~~~~~~~~~

If desired, the constituents (“members”) of each cluster can be
output to a file using ``--members=outputfilename.tab``.
Clusters are separated by empty lines and order the same as
in the clonotypes table.

In the members table, additional fields are added that are intended
to describe “mutation rates”. These fields are named
``XXX_mindiffrate``, where ``XXX`` is ``CDR3_nt``, ``CDR3_aa``, ``VDJ_nt``, and ``VDJ_aa``.

Within each cluster, the row with the lowest ``V_SHM`` value
(the least mutated V) is chosen as reference. If the ``V_SHM``
is higher than ``--v-shm-threshold``, the ``_mindiffrate`` fields are not computed.

To compute a field such as ``CDR3_nt_mindiffrate`` for a row, the edit distance between
``CDR3_nt`` of this row and of the reference row are computed and divided by
the length of ``CDR3_nt`` of the reference row (and multiplied by 100 to give a percentage).


.. _germline-filters:

Germline and pre-germline filtering
===================================

V gene sequences found by the clustering step of the program (the ``discover`` subcommand) are
stored in the ``candidates.tab`` file. The entries are “candidates” because many of these will be
PCR or other artifacts and therefore do not represent true novel V genes. The germline and
pre-germline filters take care of removing artifacts. The germline filter is the “real” filter and
used only in the last iteration in order to obtain the final gene database. The pre-germline filter
is less strict and used in all the earlier iterations.

The germline filters are implemented in the ``igdiscover germlinefilter`` subcommand. It performs the
following filtering and processing steps:

* Discard sequences with ``N`` bases
* Discard sequences that come from a consensus over too few source sequences
* Discard sequences with too few unique CDR3s (CDR3s_exact column)
* Discard sequences with too few unique Js (Js_exact column)
* Discard sequences identical to one of the database sequences (if DB given)
* Discard sequences that do not match a set of known good motifs
* Discard sequences that contain a stop codon (has_stop column)
* Discard duplicate sequences
* Discard cross-mapping artifacts
* Discard sequences whose “allele ratio” is too low.

If a whitelist of sequences is provided (by default, this is the input V gene database), then the
candidates that appear on it

* are not checked for the cluster size criterion,
* do not need to match a set of known good motifs,
* are never considered duplicates (but they are checked for
  cross-mapping and for the allele ratio),
* are allowed to contain a stop codon.

Whitelisting allows IgDiscover to identify known germline sequences that are expressed at low
levels in a library. If enabled with ``whitelist: true`` (the default) in the pregermline and
germline filter sections of the configuration file, the sequences present in the starting database
are treated as validated germline sequences and will not be discarded due to too small cluster
size as long as they fulfill the remaining criteria (unique_cdr3s, unique_js etc.).

You can see why a candidate was filtered by inspecting the
:ref:`annotated_V_*.tab files <annotated_v_tab>`


.. _cross-mapping:

Cross-mapping artifacts
-----------------------

If two very similar sequences appear in the database used by IgBLAST,
then sequencing errors may lead to one sequence incorrectly being assigned
to the other. This is particularly problematic if one of the sequences is
highly expressed while the other is not expressed at all. The not expressed
sequence is even included in the list of V gene candidates because it is
in the input database and therefore whitelisted. We call this a “cross-mapping
artifact”.

The germline filtering step of IgDiscover therefore aims to eliminate
cross-mapping artifacts by checking all pairs of sequences for the following:

* The two sequences have a distance of 1,
* they are both in the database for that particular iteration (only then
  can cross-mapping occur)
* the ratio between the expression levels of the two sequences (using
  the cluster_size field in the ``candidates.tab`` file) is less than the value
  ``cross_mapping_ratio`` defined in the configuration file (0.02 by default).

If all that is the case, then the sequence with the lower expression is
discarded.


.. _allele-ratio:

Allele-ratio filtering
----------------------

When multiple alleles of the same gene appear in the list of V gene candidates,
such as IGHV1-2*02 and IGHV1-2*04, the germline filter computes the ratio
of the values in the ``exact`` and the ``clonotypes`` columns between them.
If the ratio is under the configured threshold, the candidate with the lower
count is discarded. See the ``exact_ratio`` and ``clonotype_ratio``
settings in the ``germline_filter`` and ``pregermline_filter`` sections
of the configuration file.


.. versionadded:: 0.7.0


Data from the Sequence Read Archive (SRA)
=========================================

To work with datasets from the Sequence Read Archive, you may want to use the
tool ``fastq-dump``, which can download the reads in the format required by
IgDiscover. You just need to know the accession number, such as “SRR2905710” and
then run this command to download the files to the current directory::

    fastq-dump --split-files --gzip SRR2905710

The ``--split-files`` option ensures that the paired-end reads are stored in two
separate files, one for the forward and one for the reverse read, respectively.
(If you do not provide it, you will get an interleaved FASTQ file that currently
cannot be read by IgDiscover). The ``--gzip`` option creates compressed output.
The command creates two files in the current directory. In the above example,
they would be named ``SRR2905710_1.fastq.gz`` and ``SRR2905710_2.fastq.gz``.

The program ``fastq-dump`` is part of the SRA toolkit. On Debian-derived
Linux distributions, you can typically install it with ``sudo apt-get install
sra-toolkit``. On Conda, install it with ``conda install -c bioconda sra-tools``.


Does random subsampling influence results?
==========================================

Random subsampling indeed influences somewhat which sequences are found by the cluster analysis,
particularly in the beginning. However, the probability is large that all highly expressed
sequences are represented in the random sample. Also, due to the database growing with subsequent
iterations, the set of sequences assigned to a single database gene becomes smaller and more
homogeneous. This makes it increasingly likely that also sequences expressed at lower levels
result in a cluster since they now make up a larger fraction of each subsample.

Also, many of the clusters which are captured in one subsample but not in the other are artifacts
that are then filtered out anyway by the pre-germline or germline filter.

On human data with a nearly complete starting database, the subsampling seems to have no influence
at all, as we determined experimentally. We repeated a run of the program four
times on the same human dataset, using identical parameters each time except that the subsampling
was done in a different way. Although intermediate results differed, all four personalized
databases that the program produced were exactly identical.

Concordance is lower, though, when the input database is not as complete as the human one.

The way in which random subsampling is done is modified by the ``seed`` configuration setting,
which is set to 1 by default. If its value is the same for two different runs of the program with
otherwise identical settings, the numbers chosen by the random number generator will be the same
and therefore also subsampling will be done in an identical way. This makes runs of the program
reproducible. In order to test how results differ when subsampling is done in a different way,
change the ``seed`` to a different value.


Logging the program’s output to a file
======================================

When you report a bug or unusual behavior to us, we might ask you to send us the output of
``igdiscover run``. You can send its output to a file by running the program like this::

    igdiscover run >& logfile.txt

And here is how to send the logging output to a file *and* also see the output in your terminal
at the same time (but you lose the colors)::

  igdiscover run |& tee logfile.txt


.. _caching:

Caching of IgBLAST results and of merged reads
==============================================

Sometimes you may want to re-analyze a dataset multiple times with different filter settings.
To speed this up, IgDiscover can cache the results of two of the most time-consuming
steps, read-merging with PEAR and running IgBLAST.

The cache is disabled by default as it uses a lot of disk space. To enable the cache, create
a file named ``~/.config/igdiscover.conf`` with the following contents::

    use_cache: true

If you do so, a directory named ``~/.cache/igdiscover/`` is created the next time you run
IgDiscover and all IgBLAST results as well as merged reads from PEAR are stored there. On
subsequent runs, the existing result is used directly without calling the respective
program, which speeds up the pipeline considerably.

The cache is only used when we are certain that the results will indeed be the same. For
example, if the IgBLAST program version or th V/D/J database changes, the cached result
is not used.

The files in the cache are compressed, but the cache may still get large over time. You can
delete the cache with ``rm -r ~/.cache/igdiscover`` to free the space.

You should also delete the cache when updating to a newer IgBLAST version as the old results
will not be used anymore.


Terms
=====

Analysis directory
    The directory that was created with ``igdiscover init``. Thus, when you use
    ``igdiscover init myexperiment``, the analysis directory is ``myexperiment/``.
    Separate analysis directories need to be created for each sample.

Starting database
    The initial list of V/D/J genes. These are expected to be in FASTA format and are copied into
    the ``database/`` directory within each analysis directory.
