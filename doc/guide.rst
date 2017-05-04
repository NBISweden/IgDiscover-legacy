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
  If you are on OS X, it may be necessary to run ``export SHELL=/bin/bash`` before continuing.

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
   Even if have have only light chains in your data, a ``D.fasta`` file needs to be provided,
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
   you may :ref:`need to adjust <configuration>`. In particular, the number of discovery rounds
   is set to 3 by default, which takes a long time. Reducing this to 2 or even 1 often works just
   as well.

3. Run the analysis

   Change into the newly created analysis directory and run the analysis::

       igdiscover run

   Depending on the size of your library, your computer, and the number of iterations, this will now take from a few hours to a day.


.. _obtaining-database:

Obtaining a V/D/J database
==========================

We use the term “database” to refer to three FASTA files that contain the sequences for the V, D
and J genes.
IMGT provides `sequences for download <http://www.imgt.org/vquest/refseqh.html>`_.
For discovering new VH genes, for example, you need to get the IGHV, IGHD and IGHJ files of your species.
As IgDiscover uses this only as a starting point, using a similar species will also work.

When using an IMGT database, it is very important to change the long IMGT sequence headers to
short headers as IgBLAST does not accept the long headers. We recommend using the program
``edit_imgt_file.pl``. If you installed IgDiscover from Conda, the script is already installed and
you can run it by typing the name. It is also
`available on the IgBlast FTP site <ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/>`_.

Run it for all three downloaded files, and then rename files appropritely to make sure that they
named ``V.fasta``, ``D.fasta`` and ``J.fasta``.

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
* A random barcode (molecular identifier). This is optional. Set the
  configuration option ``barcode_length_5p`` to 0 if you don’t have random barcodes
  or if you don’t want the program to use them.
* Optionally, a run of G nucleotides. This is an artifact of the RACE protocol (Rapid
  amplification of cDNA ends). If you have this, set ``race_g`` to ``true`` in the configuration file.
* 5' UTR
* Leader
* Re-arranged V, D and J gene sequences for heavy chains; only V and J for light chains
* An optional random barcode. Set the configuration option ``barcode_length_3p`` to the length of
  this barcode. You can currently not have both a 5' and a 3' barcode.
* The reverse primer. This is optional.

We use IgBLAST to detect the location of the V, D, J genes (run as
part of the ``igdiscover igblast`` subcommand), and the remaining parts
are detected subsequently with ``igdiscover parse``. The G nucleotides
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

The main parameters parameters that may require adjusting are the following.

The ``iterations`` option sets the number of rounds of V gene discovery
that will be performed. By default, three iterations are run. Even with a very restricted
starting V database (for example with only a single V gene sequence),
this is usually sufficient to identify most novel germline sequences.

When the starting database is more complete, for example, when analyzing
a human IgM library with the current IMGT heavy chain database, a single
iteration may be sufficient to produce an individualized database.

If you do not want to discover any new genes and only want to produce an
expression profile, for example, then use ``iterations: 0``.

The ``ignore_j`` option should be set to ``true`` when producing a V gene
database for a species where J sequences are unknown::

    ignore_j: true

Setting the parameters ``stranded``, ``forward_primers`` and ``reverse_primers``
to the correct values can be used to remove 5' and 3' primers from the sequences.
Doing this is not strictly necessary for IgDiscover. It is simplest
if you do not specify any primer sequences.


Pregermline and Germline filter criteria
----------------------------------------

This provides IgDiscover with stringency requirements for V gene discovery
that enable the program to filter out false positives. Usually the ”pregermline
filter” can be used in the default mode since all these sequences will be
subsequently passed to the higher stringency ”germline filter” where the
criteria are set to maximize stringency. Here is how it looks in the configuration
file::

   pre_germline_filter:
     unique_cdr3s: 2      # Minimum number of unique CDR3s (within exact matches)
     unique_js: 2         # Minimum number of unique J genes (within exact matches)
     check_motifs: false  # Check whether 5' end starts with known motif
     whitelist: true      # Add database sequences to the whitelist
     cluster_size: 0      # Minimum number of sequences assigned to cluster
     differences: 0       # Merge sequences if they have at most this number of differences
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
     differences: 0       # Merge sequences if they have at most this number of differences
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

If the ``differences`` parameter is set to a value higher than 0, the germline filter inspects
clusters of sequences that are closely related (when the edit distance between them is at
most ``differences``) and retains only the most common sequence of each cluster. Previously, we
believed this would removes some false positives due to accumulated random sequence errors of highly
expressed alleles that otherwise would pass the cutoff criteria. However, we found out that we miss
true positives, in particular if there are two alleles in the sample that differ in only a single
nucleotide. We have now implemented other measures to avoid false positives and recommend against
setting the ``differences`` to something other than ``0``.

Read also about the :ref:`cross mapping <cross-mapping>`, for which germline filtering corrects, and
about the :ref:`germline filters <germline-filters>`.

.. versionchanged::
   The default for the ``differences`` setting was changed from 1 to 0.

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

final/database/species_(V,D,J).fasta
    These three files represent the final, individualized V/D/J database found by IgDiscover.
    The D and J files are copies of the original starting database;
    they are not updated by IgDiscover.

final/V_dendrogram.pdf
    A dendrogram of all V sequences in the individualized database.

final/igblast.txt.gz
    IgBLAST result (compressed) of running IgBLAST with the discovered database.

final/assigned.tab.gz
    V/D/J gene assignments and other information for each sequence.
    The file is created by parsing the IgBLAST output in the ``igblast.txt.gz`` file.
    This is a table that contains one row for each input sequence.
    See below for a detailed description of the columns.

final/filtered.tab.gz
    Filtered V/D/J gene assignments. This is the same as the assigned.tab file mentioned above, but with low-quality assignments filtered out.
    Run ``igdiscover filter --help`` to see the filtering criteria.

final/V_usage.tab, final/V_usage.pdf
    The V gene expression counts, derived from the IgBLAST results.
    The .tab file contains the counts as a table, while the pdf file contains a plot of the same values.

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
    This is a list of sequences found through the *windowing strategy* or *linkage cluster analysis*, as discussed in our paper.

iteration-xx/new_V_germline.fasta, iteration-xx/new_V_pregermline.fasta
    The discovered list of V genes for this iteration.
    The file is created from the ``candidates.tab`` file by applying either the germline or pre-germline filter.
    The file resulting from application of the germline filter is used in the last iteration only.
    The file resulting from application of the pre-germline filter is used in earlier iterations.


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


assigned.tab.gz
---------------

This file is a gzip-compressed table with tab-separated values.
It is created by ``igdiscover parse`` and is the result of parsing IgBLAST output.
It contains a few additional columns that do not come directly from IgBLAST.
In particular, the CDR3 sequence is detected, the sequence before the V gene match is split into *UTR* and *leader*, and
the RACE-specific run of G nucleotides is also detected.
The first row is a header row with column names.
Each subsequent row describes the IgBLAST results for a single pre-processed input sequence.

Note: This file is typically quite large.
LibreOffice can open the file directly (even though it is compressed), but make sure you have enough RAM.

Columns:

count
    How many copies of input sequence this query sequence represents. Copied from the ``;size=3;`` entry in the FASTA
    header field that is added by ``VSEARCH -derep_fulllength``.

V_gene, D_gene, J_gene
    V/D/J gene match for the query sequence

stop
    whether the sequence contains a stop codon (either “yes” or “no”)

productive

V_covered, D_covered, J_covered
    percentage of bases of the reference gene that is covered by the bases of the query sequence

V_evalue, D_evalue, J_evalue
    E-value of V/D/J hit

FR1_SHM, CDR1_SHM, FR2_SHM, CDR2_SHM, FR3_SHM, V_SHM, J_SHM
    rate of somatic hypermutation (actually, an error rate)

V_errors, J_errors
    Absolute number of errors (differences) in the V and J gene match

UTR
    Sequence of the 5' UTR (the part before the V gene match up to, but not including, the start codon)

leader
    Leader sequence (the part between UTR and the V gene match)

CDR1_nt, CDR1_aa, CDR2_nt, CDR2_aa, CDR3_nt, CDR3_aa
    nucleotide and amino acid sequence of CDR1/2/3

V_nt, V_aa
    Nucleotide and amino acid sequence of V gene match

V_CDR3_start
    Start coordinate of CDR3 within ``V_nt``. Set to zero if no CDR3 was detected.
    Comparisons involving the V gene ignore those V bases that are part of the CDR3.

V_end, VD_junction, D_region, DJ_junction, J_start
    nucleotide sequences for various match regions

name, barcode, race_G, genomic_sequence
    see the following explanation

The UTR, leader, barcode, race_G and genomic_sequence columns are filled in the following way.

1. Split the 5' end barcode from the sequence (if barcode length is zero, this will be empty), put it in the **barcode** column.
2. Remove the initial run of G bases from the remaining sequence, put that in the **race_G** column.
3. The remainder is put into the **genomic_sequence** column.
4. If there is a V gene match, take the sequence *before* it and split it up in the following way. Search for the start codon and write the part before it into the **UTR** column. Write the part starting with the start column into the **leader** column.


filtered.tab.gz
---------------

This table is the same as the ``assigned.tab.gz`` table, except that rows containing low-quality matches have been filtered out.
Rows fulfilling any of the following criteria are filtered:

- The J gene was not assigned
- A stop was codon found
- The V gene coverage is less than 90%
- The J gene coverage is less than 60%
- The V gene E-value is greater than 10\ :sup:`-3`


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
    The number of exact occurrences of the consensus sequence among all sequences assigned to the source gene.

    To clarify, we describe how the set of exact sequences is found:
    First, all sequences assigned to a source gene are clustered.
    A consensus is then computed from each cluster.
    Then we look back at *all* sequences assigned to the source gene and find exact occurrences of that consensus sequence.

Js_exact
    How many unique J genes were used by the sequences in the set of exact sequences (described above).

CDR3s_exact
    How many unique CDR3 sequences were used by the sequences in the set of exact sequences (described above).

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

approx
    Number of approximate occurrences of the consensus sequence among all sequences assigned to the source gene.
    See the description for the *exact* column.
    This *approximate set* is similar to the *exact set*, except that a difference up to a given percentage is allowed when comparing the consensus sequence to the other sequences.

Js_approx
    Same as *Js_exact*, except that it refers to the approximate occurrences of the consensus sequence.

CDR3s_approx
    Same as *CDR3s_exact*, except that it refers to the approximate occurrences of the consensus sequence.


.. _gene-names:

Novel V gene names
-------------------

Each V gene discovered by IgDiscover gets a unique name such as “VH4.11_S1234”.
The “VH4.11” is the name of the database gene to which the novel
V gene was initially assigned. The number *1234* is derived from the nucleotide
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


The following subcommands are used internally, and listed here for completeness.

parse
    Parse IgBLAST output and write out a tab-separated table

filter
    Filter table with parsed IgBLAST results

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


.. _germline-filters:

Germline and pre-germline filtering
===================================

V gene sequences found by the clustering step of the program (the ``discover`` subcommand) are
stored in the ``candidates.tab`` file. The entries are “candidates” because many of these will be
PCR or other artifacts and therefore do not represent true novel V genes. The germline and
pre-germline filters take care of removing artifacts. They germline filter is the “real” filter and
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
* Discard near-duplicate sequences
* Discard cross-mapping artifacts
* Discard sequences whose “allele ratio” is too low.

If a whitelist of sequences is provided (by default, this is the input V gene database), then the
candidates that appear on it

* are not checked for the cluster size criterion,
* do not need to match a set of known good motifs,
* are never considered near-duplicates (but they are checked for
  cross-mapping and for the allele ratio),
* are allowed to contain a stop codon.

Whitelisting allows IgDiscover to identify known germline sequences that are expressed at low
levels in a library. If enabled with ``whitelist: true`` (the default) in the pregermline and
germline filter sections of the configuration file, the sequences present in the starting database
are treated as validated germline sequences and will not be discarded if due to too small cluster
size as long as they fulfill the remaining criteria (unique_cdr3s, unique_js etc.).


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
of ``CDR3s_exact`` between them. If the ratio is under a threshold, the
lower-expressed candidate is discarded. The default threshold is 0.1 and can
be modified in the configuration file by adjusting the ``allele_ratio``
settings within the germline filter sections.


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


Terms
=====

Analysis directory
    The directory that was created with ``igdiscover init``. Separate ones are created for
    each experiment. When you used ``igdiscover init myexperiment``, the analysis directory
    would be ``myexperiment/``.

Starting database
    The initial list of V/D/J genes. These are expected to be in FASTA format and are copied into
    the ``database/`` directory within each analysis directory.
