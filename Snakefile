"""
Required modules:

module load usearch
module load cutadapt

module use -a /proj/b2013006/sw/modules
module load snakemake
module load flash
"""

USEARCH = '/proj/b2011210/dlbin/usearch7.0.1001_i86linux64'

rule merge_persson:
	output: "data/persson-merged.fastq.gz"
	input: expand("raw/lane1_Undetermined_L001_R{r}_001.fastq", r=(1,2))
	resources: time=60
	threads: 8
	shell:
		# -M: maximal overlap (2x300, 420-450bp expected fragment size)
		"flash -t {threads} -c -M 100 {input} | pigz > {output}"

rule read_length_histogram:
	output:
		txt="stats/persson-merged-readlengthhisto.txt", 
		pdf="stats/persson-merged-readlengthhisto.pdf"
	input:
		"data/persson-merged.fastq.gz"
	shell:
		"sqt-fastqhisto --plot {output.pdf} {input}  > {output.txt}"

rule demultiplex:
	# TODO
	# This isn’t doing real demultiplexing right now since the data
	# isn’t actually multiplexed.
	output:	fastq="data/persson-demultiplexed.fastq.gz"
	input: fastq="data/persson-merged.fastq.gz"
	resources: time=60
	shell:
		# This is the forward primer. The reverse primer is: GTAGTCCTTGACCAGGCAGCCCAG
		"cutadapt -g ^GCCCAGGTGAAACTGCCTCGAG --discard-untrimmed -o {output.fastq} {input.fastq}"

rule usearch_fastq_to_fasta:
	output: fasta="data/{dataset}.fasta"
	input: fastq="data/{dataset}.fastq"
	shell:
		"{USEARCH} -fastq_filter {input.fastq} -fastq_minlen 400 -fastq_maxee 1 -fastaout {output.fasta}"

rule uncompress:
	output: "{file}"
	input: "{file}.gz"
	shell: "zcat {input} > {output}"
