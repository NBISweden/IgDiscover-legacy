# kate: syntax Python;
"""
Required modules:

module use /proj/b2014342/sw/modules
module load abpipe

Each stage of the pipeline is in a different directory. They are:

1. reads/ -- raw reads
2. merged/ -- paired-end reads merged into single reads
3. trimmed/ -- adapters removed
4. filtered/ -- FASTQ converted to FASTA with low-quality reads removed
5. unique/ -- duplicate reads collapsed into one
6. clustered/ -- clustered into groups of similar sequences
7. assigned/ -- IgBLAST assignment
"""
from sqt.dna import reverse_complement

include: "pipeline.conf"

# This command is run before every shell command and helps to catch errors early
shell.prefix("set -euo pipefail;")

rule all:
	input:
		expand("fastqc/reads.{r}.zip", r=(1, 2)),
		"stats/readlengthhisto.pdf",
		"clustered.fasta",
		"igblast.txt",


#rule uncompress_reads:
	#input: "reads/{file}.fastq.gz"
	#output: "reads/{file}.fastq"
	#shell:
		#"zcat {input} > {output}"


if MERGE_PROGRAM == 'flash':
	rule flash_merge:
		"""Use FLASH to merge paired-end reads"""
		output: "merged.fastq.gz"
		input: "reads.1.fastq", "reads.2.fastq"
		resources: time=60
		threads: 8
		shell:
			# -M: maximal overlap (2x300, 420-450bp expected fragment size)
			"flash -t {threads} -c -M {FLASH_MAXIMUM_OVERLAP} {input} | pigz > {output}"
elif MERGE_PROGRAM == 'pear':
	rule pear_merge:
		"""Use pear to merge paired-end reads"""
		output: fastq="merged.fastq.gz"
		input: "reads.1.fastq", "reads.2.fastq"
		resources: time=60
		threads: 8
		shell:
			r"""
			pear -f {input[0]} -r {input[1]} -o pear && \
			gzip < pear.assembled.fastq > {output.fastq} && rm pear.assembled.fastq
			"""
else:
	sys.exit("MERGE_PROGRAM not recognized")


rule read_length_histogram:
	# TODO on which data should this be created?
	output:
		txt="stats/readlengthhisto.txt",
		pdf="stats/readlengthhisto.pdf"
	input:
		fastq="merged.fastq.gz"
	shell:
		"sqt-fastqhisto --plot {output.pdf} {input}  > {output.txt}"


# Adjust the primer sequences so they are correctly reverse-complemented.
# When a forward primer fwd and a reverse primer rev are given, then we need to
# search for:
# * fwd in the beginning, revcomp(rev) in the end
# * rev in the beginning, revcomp(fwd) in the end
PRIMERS = [
	FORWARD_PRIMERS[:], [reverse_complement(seq) for seq in REVERSE_PRIMERS]
]
PRIMERS[0].extend(REVERSE_PRIMERS)
PRIMERS[1].extend(reverse_complement(seq) for seq in FORWARD_PRIMERS)

rule trim_primers:
	output: fastq="trimmed.fastq.gz"
	input: fastq="merged.fastq.gz"
	resources: time=60
	params:
		five_p=" ".join("-g ^{}".format(seq) for seq in PRIMERS[0]),
		three_p=" ".join("-a {}$".format(seq) for seq in PRIMERS[1])
	shell:
		r"""
		cutadapt --discard-untrimmed {params.five_p} {input.fastq} | \
		cutadapt --discard-untrimmed {params.three_p} -o {output.fastq} -
		"""


rule fastqc:
	output:
		zip='fastqc/{file}.zip',
		png='fastqc/{file}/Images/per_base_quality.png',
		html='fastqc/{file}_fastqc.html'
	input: fastq='{file}.fastq'
	shell:
		r"""
		rm -rf fastqc/{wildcards.file}/ fastqc/{wildcards.file}_fastqc/ && \
		fastqc -o fastqc {input} && \
		mv fastqc/{wildcards.file}_fastqc.zip {output.zip} && \
		unzip -o -d fastqc/ {output.zip} && \
		mv fastqc/{wildcards.file}_fastqc/ fastqc/{wildcards.file}
		"""


rule fastq_to_fasta:
	"""
	* Convert from FASTQ to FASTA
	* Remove low-quality sequences
	* Discard too short sequences
	"""
	output: fasta="filtered.fasta"
	input: fastq="trimmed.fastq.gz"
	params: max_errors="--max-errors {MAXIMUM_EXPECTED_ERRORS}" if MAXIMUM_EXPECTED_ERRORS is not None else ""
	shell:
		"sqt-fastqmod {params.max_errors} --minimum-length {MINIMUM_MERGED_READ_LENGTH} --fasta {input.fastq} > {output.fasta}"


rule dereplicate:
	"""Dereplication with usearch"""
	output: fasta="unique.fasta"
	input: fasta="filtered.fasta"
	shell:
		"""vsearch --derep_fulllength {input.fasta} --strand both --output {output.fasta} --sizeout"""


rule cluster:
	"""
	TODO Daniel ran this three times (at 99%, 98% and 97% identity) in order to
	avoid a specific type of misclustering.
	"""
	output:
		fasta="clustered.fasta",  # centroids
		uc="clustered.uc"
	input: fasta="unique.fasta"
	resources: time=36*60, mem=32000
	threads: 4
	shell:
		# TODO -idprefix 5?
		r"""
		{CLUSTER_PROGRAM} -threads {threads} -cluster_fast {input.fasta} -id 0.97 -uc {output.uc} \
			-idprefix 5 -sizeout --centroids {output.fasta}
		"""


rule makeblastdb:
	input: fasta="database/{organism}_{gene}.fasta"
	output: "database/{organism}_{gene}.nhr" # and nin nog nsd nsi nsq
	params: dbname="database/{organism}_{gene}"
	shell:
		r"""
		makeblastdb -parse_seqids -dbtype nucl -in {input.fasta} -out {params.dbname}
		"""


rule igblast:
	output:
		txt="igblast.txt"
	input:
		fasta="unique.fasta",
		db_v="database/{species}_V.nhr".format(species=SPECIES),
		db_d="database/{species}_D.nhr".format(species=SPECIES),
		db_j="database/{species}_J.nhr".format(species=SPECIES)
	shell:
		#-auxiliary_data $IGDATA/optional_file/{SPECIES}_gl.aux
		r"""
		igblastn \
			-germline_db_V database/{SPECIES}_V \
			-germline_db_D database/{SPECIES}_D \
			-germline_db_J database/{SPECIES}_J \
			-organism {SPECIES} \
			-ig_seqtype Ig \
			-num_threads 1 \
			-domain_system imgt \
			-num_alignments_V 1 \
			-num_alignments_D 1 \
			-num_alignments_J 1 \
			-out {output.txt} \
			-query {input.fasta} \
			-outfmt '7 qseqid qstart qseq sstart sseq pident'
		"""


rule igblastwrp:
	"""TODO run it on clustered sequences, not on dereplicated ones (?)"""
	output: table="igblastwrp-table.L2.txt"
	input: "filtered.fasta"
	resources: time=60
	threads: 8
	shell:
		r"""
		igblastwrp -R {RECEPTOR_CHAIN} -S {SPECIES} -p {threads} {input} igblast-table
		"""


rule ungzip:
	output: "{file}.fastq"
	input: "{file}.fastq.gz"
	shell: "zcat {input} > {output}"
