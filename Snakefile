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

include: "pipeline.conf"


# This command is run before every shell command and helps to catch errors early
shell.prefix("set -euo pipefail;")

rule all:
	input:
		expand("fastqc/reads.{r}.zip", r=(1, 2)),
		"stats/readlengthhisto.pdf",
		"clustered.fasta",
		"igblast-table.L2.txt",


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


rule trim_primers:
	output: fastq="trimmed.fastq.gz"
	input: fastq="merged.fastq.gz"
	resources: time=60
	params:
		five_p=" ".join("-g ^{}".format(seq) for seq in FIVE_PRIME_PRIMERS),
		three_p=" ".join("-a {}".format(seq) for seq in THREE_PRIME_PRIMERS)
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


rule usearch_fastq_to_fasta:
	"""
	Convert from FASTQ to FASTA; remove low-quality sequences;
	discard too short sequences.
	"""
	output: fasta="filtered.fasta"
	input: fastq="trimmed.fastq"
	shell:
		"{USEARCH} -fastq_filter {input.fastq} -fastq_minlen 400 -fastq_maxee {MAXIMUM_EXPECTED_ERRORS} -fastaout {output.fasta}"
		# alternatively, this should work:
		# "sqt-fastqmod --max-errors {MAXIMUM_EXPECTED_ERRORS} --minimum-length 400 --fasta {input.fastq} > {output.fasta}"


rule usearch_derep_fulllength:
	"""Dereplication with usearch"""
	output: fasta="unique.fasta"
	input: fasta="filtered.fasta"
	shell:
		"""{USEARCH} -derep_fulllength {input.fasta} -output {output.fasta} -sizeout"""


rule usearch_cluster:
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
		# TODO uses only 1 thread
		r"""
		{USEARCH} -threads {threads} -cluster_fast {input.fasta} -id 0.97 -uc {output.uc} \
			-idprefix 5 --centroids {output.fasta} -sizeout
		"""


rule igblast:
	"""TODO run it on clustered sequences, not on dereplicated ones (?)"""
	output: table="igblast-table.L2.txt", igblast="igblast.txt"
	input: "filtered.fasta"
	resources: time=60
	threads: 8
	shell:
		r"""
		igblastwrp --debug -R {RECEPTOR_CHAIN} -S {SPECIES} -p {threads} {input} igblast-table > {output.igblast}
		"""


rule ungzip:
	output: "{file}.fastq"
	input: "{file}.fastq.gz"
	shell: "zcat {input} > {output}"
