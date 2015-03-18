# kate: syntax Python;
"""
Required modules:

module use /proj/b2014342/sw/modules
module load abpipe

These are the main files that the pipeline creates, in the order in which they
are created:

merged.fastq.gz -- merged reads
trimmed.fastq.gz -- primers removed from merged reads
filtered.fasta  -- too short sequences removed, converted to FASTA
unique.fasta -- collapsed sequences
igblast.txt -- raw IgBLAST output
table.tab -- result of parsing IgBLAST output
"""
from sqt.dna import reverse_complement

# Set defaults for some configurable values.

# Limit processing to this number of reads. Set to None to process all.
LIMIT = None

# Which program to use for clustering. This can be either a path to the
# usearch binary or vsearch.
CLUSTER_PROGRAM = 'vsearch'

# Filter out reads that have more than this number of expected errors.
# Set to None to disable.
MAXIMUM_EXPECTED_ERRORS = None

# Whether to trim primers. Can be set to True or False.
TRIM_PRIMERS = True

try:
	include: "pipeline.conf"
except WorkflowError:
	sys.exit("Pipeline configuration file 'pipeline.conf' not found. Please create it!")

if not REVERSE_PRIMERS:
	sys.exit("The list of REVERSE_PRIMERS is empty. This will currently not work correctly.")

# This command is run before every shell command and helps to catch errors early
shell.prefix("set -euo pipefail;")

rule all:
	input:
		#expand("fastqc/reads.{r}.zip", r=(1, 2)),
		"stats/readlengthhisto.pdf",
		"clustered.fasta",
		"table.tab",
		"v_usage.tab",
		"v_usage.pdf"


if MERGE_PROGRAM == 'flash':
	rule flash_merge:
		"""Use FLASH to merge paired-end reads"""
		output: "merged.fastq.gz"
		input: "reads.1.fastq", "reads.2.fastq"
		resources: time=60
		threads: 8
		log: 'flash.log'
		shell:
			# -M: maximal overlap (2x300, 420-450bp expected fragment size)
			"flash -t {threads} -c -M {FLASH_MAXIMUM_OVERLAP} {input} 2> >(tee {log} >&2) | pigz > {output}"
elif MERGE_PROGRAM == 'pear':
	rule pear_merge:
		"""Use pear to merge paired-end reads"""
		output: fastq="merged.fastq.gz"
		input: "reads.1.fastq", "reads.2.fastq"
		resources: time=60
		threads: 8
		shell:
			r"""
			pear -j {threads} -f {input[0]} -r {input[1]} -o pear && \
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

rule trim_primers_fivep:
	output: fastq="trimmed-fivep.fastq.gz"
	input: fastq="merged.fastq.gz"
	resources: time=60
	log: 'cutadapt-fiveprime.log'
	params:
		five_p=" ".join("-g ^{}".format(seq) for seq in PRIMERS[0])
	shell:
		r"""
		cutadapt --discard-untrimmed {params.five_p} -o {output.fastq} {input.fastq} | tee {log}
		"""


rule trim_primers_threep:
	output: fastq="trimmed.fastq.gz"
	input: fastq="trimmed-fivep.fastq.gz"
	resources: time=60
	log: 'cutadapt-threeprime.log'
	params:
		three_p=" ".join("-a {}$".format(seq) for seq in PRIMERS[1])
	shell:
		r"""
		cutadapt --discard-untrimmed {params.three_p} -o {output.fastq} {input.fastq} | tee {log}
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
	input: fastq="trimmed.fastq.gz" if TRIM_PRIMERS else "merged.fastq.gz"
	params: max_errors="--max-errors {MAXIMUM_EXPECTED_ERRORS}" if MAXIMUM_EXPECTED_ERRORS is not None else ""
	shell:
		"sqt-fastqmod {params.max_errors} --minimum-length {MINIMUM_MERGED_READ_LENGTH} --fasta {input.fastq} > {output.fasta}"


rule dereplicate:
	"""Collapse identical sequences with VSEARCH"""
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
	log: 'database/{organism}_{gene}.log'
	threads: 100  # force to run as only job
	shell:
		r"""
		makeblastdb -parse_seqids -dbtype nucl -in {input.fasta} -out {params.dbname} >& {log}
		grep '^Error:' {log} && {{ echo "makeblastdb failed when creating {params.dbname}"; false; }} || true
		"""


rule abpipe_igblast:
	output:
		txt="igblast.txt"
	input:
		fasta="unique.fasta",
		db_v="database/{species}_V.nhr".format(species=SPECIES),
		db_d="database/{species}_D.nhr".format(species=SPECIES),
		db_j="database/{species}_J.nhr".format(species=SPECIES)
	params:
		limit='--limit {}'.format(LIMIT) if LIMIT else ''
	threads: 16
	shell:
		#-auxiliary_data $IGDATA/optional_file/{SPECIES}_gl.aux
		r"""
		abpipe igblast --threads {threads} {params.limit} --species {SPECIES} database/ {input.fasta} > {output.txt}
		"""


rule parse_igblast:
	output:
		tab="table.tab"
	input:
		txt="igblast.txt",
		fasta="unique.fasta"
	params:
		dirname=os.path.basename(os.getcwd())
	shell:
		"abpipe parse --rename {params.dirname}_ {input.txt} {input.fasta} > {output.tab}"


rule count_and_plot:
	output:
		plot="v_usage.pdf",
		counts="v_usage.tab"
	input:
		tab="table.tab"
	shell:
		"abpipe count {input.tab} {output.plot} > {output.counts}"


rule ungzip:
	output: "{file}.fastq"
	input: "{file}.fastq.gz"
	shell: "zcat {input} > {output}"
