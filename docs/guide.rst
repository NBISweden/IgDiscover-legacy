==========
User guide
==========


IgDiscover overview
===================

IgDiscover works on a single library at a time. It creates a subdirectory for
the library, which contains all intermediate and result files.

To start an analysis, you need:

1. Two FASTQ files with paired-end reads
2. A database of V/D/J genes (three FASTA files named ``v.fasta``, ``d.fasta``, ``j.fasta``)
3. A configuration file that describes the library

To run an analysis, proceed as follows.

1. Create and initialize the analysis directory.

   First, pick a name for your analysis. We will use ``myexperiment`` in the following.
   Run ``igdiscover init``::

       igdiscover init myexperiment

   A dialog will appear and ask for the file with the *first* reads.
   Find the ``mylibrary.1.fastq.gz`` file or whatever it is called and select it.
   You do not need to choose the second read file!
   It is found automatically.

   Next, choose the directory with your database.
   The directory must contain the three files ``v.fasta``, ``d.fasta``, ``j.fasta``.
   These files contain the VH, DH, JH gene sequences, respectively.

   If you do not want a graphical user interface, use the two command-line
   parameters ``--db`` and ``--reads`` to provide this information instead::

       igdiscover init --db path/to/my/database/ --reads mylibrary.1.fastq.gz myexperiment

   In any case, an analysis directory named ``myexperiment`` will have been created.

2. Adjust the configuration file

   The previous step created a configuration file named ``myexperiment/igdiscover.yaml``.
   Use a text editor to modify the file.
   The settings you should make sure to check are:

   - ``iterations`` and
   - ``species``

   The ``iterations`` specifies the number of rounds of V gene discovery that will be performed.
   Four or five rounds are usually sufficient.
   You can also set this to zero, in which case the initial database will be used unchanged.

   The ``species`` must be set to a value supported by IgBLAST.
   As IgDiscover uses custom V/D/J databases, this setting does not influence results too much.

   Setting the parameters ``stranded``, ``forward_primers`` and ``reverse_primers`` to the correct values can be used to remove 5' and 3' primers from the sequences.
   Doing this is not strictly necessary for IgDiscover!
   It is simplest if you do not specify any primer sequences.

3. Run the analysis

   Change into the newly created analysis directory and run the analysis::

       igdiscover run

   This will use `snakemake <http://snakemake.bitbucket.org/>`_ to run the pipeline on all available CPU cores.

   Add the ``-p`` parameter to print out slightly more information::

       igdiscover run -p

   Depending on the size of your library, your computer, and the number of iterations, this will now take from a few hours to a day.



The analysis directory
======================

IgDiscover writes all intermediate files, the final V gene database, statistics and plots into the analysis directory that was created with ``igdiscover init``.
The files in the ``final/`` subdirectory are likely the most relevant ones.

These are the files and subdirectories can be found in the directory.
(Subdirectories are described in detail below.)

igdiscover.yaml
    The configuration file.
    Make sure to adjust this to your needs as described above.

reads.1.fastq.gz
reads.2.fastq.gz
    Symbolic links to the raw, paired-end input reads.

database/
    The input V/D/J database (simply three FASTA files).
    The files are a copy of the ones you selected when running ``igdiscover init``.

reads/
    Processed reads (merged etc.)

iteration-xx/
    Iteration-specific analysis directory.
    Each iteration is run in one of these directories.
    ``iteration-01`` contains the first iteration, ``iteration-02`` the second etc.
    Iteration 01 uses the original input database (the one also found in the ``database/`` directory),
    and creates an updated database from it (with new V genes).
    The new database is used as input for the next iteration.

final/
    After the last iteration, IgBLAST is run again on the input sequences, but using the final database (the one created in the very last iteration).
    This directory contains all the results, such as plots of the repertoire profiles.
    If you set the number of iterations to 0 in the configuration file, this directory is the only one that is created.


Final results
-------------

Final results can be found in the ``final`` subdirectory of the analysis directory.

final/database/species_V.fasta
final/database/species_D.fasta
final/database/species_J.fasta
    The final, individualized V, D, J database found by IgDiscover.
    The D and J files are copies of the original starting database;
    they are not updated by IgDiscover.

final/V_dendrogram.pdf
    A dendrogram of all V sequences in the individualized database.

final/unique.igblast.txt.gz
    IgBLAST result (compressed) of running IgBLAST with the discovered database.

final/unique.assigned.tab.gz
    Parsed IgBLAST results as a table (created from the igblast.txt.gz file).
    This table contains one row for each input sequence.
    See below for a detailed description of the columns.

final/unique.filtered.tab.gz

final/V_usage.tab
final/V_usage.pdf
    The V gene expression counts, derived from the IgBLAST results.
    The .tab file contains the counts as a table, while the pdf file contains a plot of the same values.

final/unique.errorhistograms.pdf

final/clusterplots/
    VH7.21_S4259.png


If you are interested in the results of each iteration, you can inspect the iteration-xx/ directories.
They are structured in the same way as the final/ subdirectory, except that the results are based on the intermediate databases of that iteration.
They also contain the following additional files.

iteration-xx/candidates.tab
iteration-xx/new_V_database.fasta



reads/merged.fastq.gz
    Reads merged with PEAR or FLASH
reads/trimmed.fastq.gz
    Merged reads with 5' and 3' primer sequences removed.
reads/filtered.fasta
    Merged, primer-trimmed sequences converted to FASTA, and too short sequences removed.
reads/unique.fasta
    Filtered sequences without duplicates (using VSEARCH)


stats/unique.readlengths.txt
stats/merged.readlengths.pdf
stats/unique.readlengths.pdf
stats/merged.readlengths.txt
stats/barcodes.txt


Creating a new IgBLAST database
===============================

To reduce confusion, do not modify the database, but create a new version
instead. Create a new directory in ``igdiscover/databases``.
Then copy FASTA files with V, D, J sequences into the directory. The files need
to be named ``rhesus_monkey_V.fasta``, ``rhesus_monkey_D.fasta`` and
``rhesus_monkey_J.fasta``. The ``makeblastdb`` program will be run automatically
by the pipeline next time it runs.




Structure of each sequence
==========================

IgDiscover assumes that its input data are overlapping paired-end reads. After
merging, they should have this structure (from 5' to 3'):

* A random barcode (molecular identifier). This is optional. Set the
  configuration option ``barcode_length`` to 0 if you don’t have random barcodes
  or if you don’t want the program to use them.
* A run of G nucleotides. This is an artifact of the RACE protocol (Rapid
  amplification of cDNA ends).
* 5' UTR
* Leader
* Re-arranged V, D and J gene sequences (in that order)

We use IgBLAST to detect the location of the V, D, J genes (run as part of the
``igdiscover igblast`` subcommand), and the remaining parts are detected
subsequently with ``igdiscover parse``. The G nucleotides after the barcode are
always split off, even if no RACE protocol was used. (This should not be a
problem in practice.) The leader sequence is detected by looking for a start
codon near 60 bp upstream of the start of the V gene match.


Novel VH gene names
===================

Each novel VH gene discovered by IgDiscover gets a unique name such as
“VH4.11_S1234”. The “VH4.11” is the name of the database gene to which the novel
VH gene was initially assigned. The number *1234* is derived from the base
sequence of the novel gene. That is, if you discover the same sequence in two
different runs of the IgDiscover, or just in different iterations, the number will
be the same. This may help when manually inspecting results.

Be aware that you still need to check the sequence itself since even different
sequences can sometimes lead to the same number (a “hash collision”).


Subcommands
===========

commonv             Find common V genes between two different antibody libraries.
igblast             Run IgBLAST.
parse               Parse IgBLAST output and write out a tab-separated table.
filter              Filter table with parsed IgBLAST results
count               Count and plot V, D, J gene usage.
group               Group sequences by barcode and V/J assignment and print each group’s consensus
multidiscover       Find V gene sister sequences shared by multiple libraries.
compose             Create new V gene database from V gene candidates.
discover            Discover candidate new V genes within a single antibody library.
init                Create and initialize a new pipeline directory.
clusterplot         For each V gene, plot a clustermap of the sequences assigned to it.
errorplot           Plot histograms of differences to reference V gene
upstream            Cluster upstream sequences (UTR and leader) for each gene
dendrogram          Draw a dendrogram of sequences in a FASTA file.
rename              Rename sequences in a target FASTA file using a template FASTA file
union               Compute union of sequences in multiple FASTA files


The assigned.tab table
======================

This file is created by
The file created by ``igdiscover parse`` is written to a file named ``...assigned.tab``. It contains the results of parsing IgBLAST output. Each row describes the result for a single query sequence. The first row is a header row.

Columns
-------

* count: How many copies of input sequence this query sequence represents. Copied from the ``;size=3;`` entry in the FASTA header field that is added by ``VSEARCH -derep_fulllength``.
* V_gene, D_gene, J_gene: V/D/J gene match for the query sequence
* stop (yes/no): whether the sequence contains a stop codon
* productive
* UTR
* leader
* V_covered, D_covered, J_covered: percentage of bases of the reference gene that is covered by the bases of the query sequence
* V_evalue, D_evalue, J_evalue: E-value of V/D/J hit
* FR1_SHM, CDR1_SHM, FR2_SHM, CDR2_SHM, FR3_SHM, V_SHM, J_SHM: rate of somatic hypermutation (actually, an error rate)
* CDR1_nt, CDR1_aa, CDR2_nt, CDR2_aa, CDR3_nt, CDR3_aa
* V_nt, V_aa: nucleotide and amino acid sequence of V gene match
* V_end, VD_junction, D_region, DJ_junction, J_start: nucleotide sequences
* name
* barcode
* race_G
* genomic_sequence

The UTR, leader, barcode, race_G and genomic_sequence columns are filled in the following way.

1. Split 5' end barcode from the sequence (if barcode length is zero, this will be empty), put it in the **barcode** column.
2. Remove the initial run of G bases from the remaining sequence, put that in the **race_G** column.
3. The remainder is put into the **genomic_sequence** column.
4. If there is a V gene match, take the sequence *before* it and split it up in the following way. Search for the start codon and write the part before it into the **UTR** column. Write the part starting with the start column into the **leader** column.


The discover.tab table
======================

The output table generated by ``igdiscover discover``, named ``...discover.tab``, has the following columns:

* ``gene``: name of the V gene

Then there are sequence counts for four subsets of all the V sequences assigned to that gene. The groups are:

* ``total``: *all* sequences assigned to this gene
* ``window`: Those sequences assigned to this gene that are within the specified error rate window (with command-line options ``--left`` and ``--right``). The consensus sequence is computed from these.
* ``exact``: Those sequences assigned to this gene that are identical to the consensus (exact matches)
* ``approx``: Those sequences assigned to this gene that match the consensus approximately (the allowed error rate is by default 1%, but can be changed with the ``--error-rate`` parameter)

For each of those four groups, the following numbers are given:

* ``_seqs``: Number of sequences in this group
* ``_unique_J``: Number of unique J genes used in this group
* ``_unique_CDR3``: Number of unique CDR3 sequences used in this group

Finally, the last columns are:

* ``N_bases``: Number of `N` bases in the consensus
* ``database_diff``: Number of differences between consensus and database sequence. This is only assigned when the database sequences were provided with ``--database``.
* ``consensus``: The consensus sequence itself


Configuration
=============

forward_primers, reverse_primers: If any primer sequences are given here, then
reads that do not have the primer sequence will be discarded.

If you use an unstranded protocol, set the ``stranded`` setting to ``false``.
The pipeline will then search also reverse-complemented reads for primers.
