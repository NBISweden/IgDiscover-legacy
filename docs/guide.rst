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

   Change into the newly created analysis directory and run the ``snakemake`` tool in order to run the analysis.
   Make sure you add the ``-j`` parameter in order to use all available cores of your machine::

       snakemake -j

   The ``-p`` parameter prints out a bit more information while IgDiscover works if you would like to see it::

       snakemake -p -j

   Depending on the size of your library, your computer, and the number of iterations, this will now take from a few hours to a day.



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

filtered.fasta
flash.log
forward-primer-trimmed.cutadapt.log
forward-primer-trimmed.fastq.gz
merged.fastq.gz
pear.assembled.fastq
pear.discarded.fastq
pear.log
pear.unassembled.forward.fastq
pear.unassembled.reverse.fastq
trimmed.cutadapt.log
trimmed.fastq.gz
unique.fasta



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
