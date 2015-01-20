# kate: syntax Python;
"""
Required modules:

module purge
module load gcc bioinfo-tools
module use /proj/b2013006/sw/modules  # Load this first to avoid the broken Uppmax cutadapt and python modules.
module load FastQC cutadapt snakemake flash pear igblastwrp

# module load usearch  # unused, we need the 64 bit version

To run this, create symlinks in the reads/ directory that point to your raw data.

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
		expand("fastqc/{dataset}.{r}.zip", r=(1, 2), dataset=DATASETS),
		expand("stats/{dataset}.readlengthhisto.pdf", dataset=DATASETS),
		expand("clustered/{dataset}.fasta", dataset=DATASETS),
		expand("igblast/{dataset}.L0.txt", dataset=DATASETS),


rule create_persson_small:
	"Create a subset of the 20M Persson reads"
	output: fastq="reads/small.{r}.fastq"
	input: fastq="reads/persson.{r}.fastq"
	shell:
		"head -n 44000000 {input} | tail -n 4000000 > {output}"


#rule uncompress_reads:
	#input: "reads/{file}.fastq.gz"
	#output: "reads/{file}.fastq"
	#shell:
		#"zcat {input} > {output}"


if MERGE_PROGRAM == 'flash':
	rule flash_merge:
		"""Use FLASH to merge paired-end reads"""
		output: "merged/{dataset}.fastq.gz"
		input: expand("reads/{{dataset}}.{r}.fastq", r=(1,2))
		resources: time=60
		threads: 8
		shell:
			# -M: maximal overlap (2x300, 420-450bp expected fragment size)
			"flash -t {threads} -c -M {FLASH_MAXIMUM_OVERLAP} {input} | pigz > {output}"
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


rule trim_primers:
	output: fastq="trimmed/{dataset}.fastq.gz"
	input: fastq="merged/{dataset}.fastq.gz"
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
	input: fastq="trimmed/{dataset}.fastq"
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
	"""TODO run it on clustered sequences, not on dereplicated ones"""
	output: "igblast/{dataset}.L0.txt"
	input: "unique/{dataset}.fasta"
	#resources: time=??
	threads: 8
	shell:
		"igblastwrp -R {RECEPTOR_CHAIN} -S {SPECIES} -p {threads} {input} igblast/{wildcards.dataset}"


rule ungzip:
	output: "{file}.fastq"
	input: "{file}.fastq.gz"
	shell: "zcat {input} > {output}"
