==========
User guide
==========


Running the pipeline
====================

Each library will be processed in a separate subdirectory of the pipeline
directory, such as ``pipeline/xyz`` for a library named ``xyz``. Use the command
``igdiscover init pipeline/xyz`` to initialize the directory for a new library. When
you run it, you will be asked to select the file that contains the reads that
you want to analyze. Only select the first file; the second file will be
determined automatically.

You will then be asked for the IgBLAST database directory to use. Next, you
should adjust the configuration file ``pipeline.conf`` if necessary.

Finally, make sure you are in the ``xyz`` directory, and then run::

    snakemake -j

The ``-j`` option makes sure that as many processors as available on your machine
are used in parallel.


Creating a new IgBLAST database
===============================

To reduce confusion, do not modify the database, but create a new version
instead. Create a new directory in ``igdiscover/databases``.
Then copy FASTA files with V, D, J sequences into the directory. The files need
to be named ``rhesus_monkey_V.fasta``, ``rhesus_monkey_D.fasta`` and
``rhesus_monkey_J.fasta``. The ``makeblastdb`` program will be run automatically
by the pipeline next time it runs.


Files
=====

reads/merged.fastq.gz -- merged reads
reads/trimmed.fastq.gz -- primers removed from merged reads
reads/filtered.fasta  -- too short sequences removed, converted to FASTA
reads/unique.fasta -- collapsed sequences (duplicates removed)
unique.igblast.txt.gz -- IgBLAST output
unique.assigned.tab -- parsed IgBLAST output as a tab-delimited table
unique.filtered.tab -- filtered version of the above
groups.tab -- sequences grouped by barcode
consensus.fasta -- contains one consensus sequence for each group
consensus.igblast.txt -- consensus sequences sent through IgBLAST
consensus.assigned.tab -- parsed IgBLAST output


reads.1.fastq.gz
reads.2.fastq.gz
igdiscover.yaml
Snakefile

reads/decompressed.1.fastq
reads/decompressed.2.fastq
reads/merged.fastq.gz
reads/forward-primer-trimmed.fastq.gz
reads/trimmed.fastq.gz
reads/unique.fasta
reads/pear.unassembled.reverse.fastq
reads/pear.unassembled.forward.fastq
reads/pear.discarded.fastq
reads/pear.log
reads/filtered.fasta

iteration-05/consensus.fasta
iteration-05/unique.fasta
iteration-05/consensus.assigned.tab.gz
iteration-05/unique.correlationVJ.pdf
iteration-05/clusterplots/VH7.21_S4259.png
iteration-05/clusterplots/done
iteration-05/consensus.igblast.txt.gz
iteration-05/candidates.tab
iteration-05/new_V_database.fasta
iteration-05/unique.errorhistograms.pdf
iteration-05/counts.txt
iteration-05/stats
iteration-05/stats/groupsizes.pdf
iteration-05/groups.tab
iteration-05/unique.V_usage.pdf
iteration-05/database
iteration-05/database/rhesus_monkey_J.fasta
iteration-05/database/rhesus_monkey_D.fasta
iteration-05/database/rhesus_monkey_V.fasta  (+ .nsi/nin/snq/nhr/nog)
iteration-05/unique.V_usage.tab
iteration-05/unique.assigned.tab.gz
iteration-05/unique.filtered.tab.gz
iteration-05/unique.igblast.txt.gz
iteration-05/V_dendrogram.pdf
iteration-05/unique.consensus.log

final/consensus.V_usage.tab
final/consensus.fasta
final/unique.fasta
final/consensus.assigned.tab.gz
final/unique.correlationVJ.pdf
final/clusterplots/VH7.21_S4259.png
final/clusterplots/done
final/consensus.igblast.txt.gz
final/unique.errorhistograms.pdf
final/consensus.correlationVJ.pdf
final/counts.txt
final/stats
final/stats/groupsizes.pdf
final/groups.tab
final/unique.V_usage.pdf
final/database/rhesus_monkey_J.fasta
final/database/rhesus_monkey_D.fasta
final/database/rhesus_monkey_V.fasta
final/unique.V_usage.tab
final/consensus.filtered.tab.gz
final/unique.assigned.tab.gz
final/unique.filtered.tab.gz
final/unique.igblast.txt.gz
final/V_dendrogram.pdf
final/consensus.V_usage.pdf
final/unique.consensus.log

stats
stats/unique.readlengths.txt
stats/merged.readlengths.pdf
stats/unique.readlengths.pdf
stats/merged.readlengths.txt
stats/barcodes.txt
database
database/rhesus_monkey_J.fasta
database/rhesus_monkey_D.fasta
database/rhesus_monkey_V.fasta





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

    commonv             Find common V genes between two different antibody
                        libraries.
    igblast             Run IgBLAST.
    parse               Parse IgBLAST output and write out a tab-separated
                        table.
    filter              Filter table with parsed IgBLAST results
    count               Count and plot V, D, J gene usage.
    group               Group sequences by barcode and V/J assignment and
                        print each group’s consensus
    multidiscover       Find V gene sister sequences shared by multiple
                        libraries.
    compose             Create new V gene database from V gene candidates.
    discover            Discover candidate new V genes within a single
                        antibody library.
    init                Create and initialize a new pipeline directory.
    clusterplot         For each V gene, plot a clustermap of the sequences
                        assigned to it.
    errorplot           Plot histograms of differences to reference V gene
    upstream            Cluster upstream sequences (UTR and leader) for each
                        gene
    dendrogram          Draw a dendrogram of sequences in a FASTA file.
    rename              Rename sequences in a target FASTA file using a
                        template FASTA file
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
