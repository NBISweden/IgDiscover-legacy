# kate: syntax Python;
"""
Required modules:

module purge
module load gcc
module load bioinfo-tools
# module load usearch  # unused, we need the 64 bit version
module load cutadapt
module load FastQC

module use /proj/b2013006/sw/modules
module load snakemake
module load flash

To run this, create symlinks in the reads/ directory that point to your raw data.

Each stage of the pipeline is in a different directory. They are:

1. reads/ -- raw reads
2. merged/ -- paired-end reads merged into single reads
3. demultiplexed/ -- demultiplexed and without adapters
4. filtered/ -- FASTQ converted to FASTA with low-quality reads removed
5. unique/ -- duplicate reads collapsed into one
6. clustered/ -- clustered into groups of similar sequences
"""

MERGE_PROGRAM = 'pear'  # either 'pear' or 'flash'

FORWARD_PRIMER = 'GCCCAGGTGAAACTGCCTCGAG'
# The reverse primer is: GTAGTCCTTGACCAGGCAGCCCAG
USEARCH = '/proj/b2011210/dlbin/usearch7.0.1001_i86linux64'

# Other configuration options:
# - species
# - which chain (heavy/light) to analyze

DATASETS = [ 'persson', 'small' ]

# This command is run before every shell command and helps to catch errors early
shell.prefix("set -euo pipefail;")

rule all:
	input:
		expand("fastqc/{dataset}.{r}.zip", r=(1, 2), dataset=DATASETS),
		expand("stats/{dataset}.readlengthhisto.pdf", dataset=DATASETS),
		expand("clustered/{dataset}.fasta", dataset=DATASETS),


rule create_persson_small:
	"Create a subset of the 20M Persson reads"
	output: fastq="reads/small.{r}.fastq"
	input: fastq="reads/persson.{r}.fastq"
	shell:
		"head -n 44000000 {input} | tail -n 4000000 > {output}"


if MERGE_PROGRAM == 'flash':
	rule flash_merge:
		"""Use FLASH to merge paired-end reads"""
		output: "merged/{dataset}.fastq.gz"
		input: expand("reads/{{dataset}}.{r}.fastq", r=(1,2))
		resources: time=60
		threads: 8
		shell:
			# -M: maximal overlap (2x300, 420-450bp expected fragment size)
			"flash -t {threads} -c -M 100 {input} | pigz > {output}"
elif MERGE_PROGRAM == 'pear':
	rule pear_merge:
		"""Use pear to merge paired-end reads"""
		output: fastq="merged/{dataset}.fastq"
		input: expand("reads/{{dataset}}.{r}.fastq", r=(1,2))
		resources: time=60
		threads: 8
		shell:
			r"""
			pear -f {input[0]} -r {input[1]} -o {wildcards.dataset} && \
			mv {wildcards.dataset}.assembled.fastq {output.fastq}
			"""
else:
	sys.exit("MERGE_PROGRAM not recognized")


rule read_length_histogram:
	# TODO on which data should this be created?
	output:
		txt="stats/{dataset}.readlengthhisto.txt",
		pdf="stats/{dataset}.readlengthhisto.pdf"
	input:
		fastq="merged/{dataset}.fastq.gz"
	shell:
		"sqt-fastqhisto --plot {output.pdf} {input}  > {output.txt}"


rule demultiplex:
	# TODO
	# This isn’t doing real demultiplexing right now since the data
	# isn’t actually multiplexed.
	output:	fastq="demultiplexed/{dataset}.fastq.gz"
	input: fastq="merged/{dataset}.fastq.gz"
	resources: time=60
	shell:
		"cutadapt -g ^{FORWARD_PRIMER} --discard-untrimmed -o {output.fastq} {input.fastq}"


rule fastqc:
	output:
		zip='fastqc/{file}.zip',
		png='fastqc/{file}/Images/per_base_quality.png',
		html='fastqc/{file}_fastqc.html'
	input: fastq='reads/{file}.fastq'
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
	output: fasta="filtered/{dataset}.fasta"
	input: fastq="demultiplexed/{dataset}.fastq"
	shell:
		"{USEARCH} -fastq_filter {input.fastq} -fastq_minlen 400 -fastq_maxee 1 -fastaout {output.fasta}"


rule usearch_derep_fulllength:
	"""Dereplication with usearch"""
	output: fasta="unique/{dataset}.fasta"
	input: fasta="filtered/{dataset}.fasta"
	shell:
		"""{USEARCH} -derep_fulllength {input.fasta} -output {output.fasta} -sizeout"""


rule usearch_cluster:
	"""
	TODO Daniel ran this three times (at 99%, 98% and 97% identity) in order to
	avoid a specific type of misclustering.
	"""
	output:
		fasta="clustered/{dataset}.fasta",  # centroids
		uc="clustered/{dataset}.uc"
	input: fasta="unique/{dataset}.fasta"
	resources: time=12*60
	#threads: 1
	shell:
		# TODO -idprefix 5?
		# TODO even with -threads {threads}, it uses only 1 thread
		r"""
		{USEARCH} -cluster_fast {input.fasta} -id 0.97 -uc {output.uc} \
			-idprefix 5 --centroids {output.fasta} -sizeout
		"""

rule ungzip:
	output: "{file}.fastq"
	input: "{file}.fastq.gz"
	shell: "zcat {input} > {output}"
