#!/bin/bash
# Run this within an activated igdiscover environment
set -euo pipefail
set -x
unset DISPLAY

pytest

rm -rf testrun
mkdir testrun
[[ -L testdata ]] || ln -s igdiscover-testdata testdata


# Test whether specifying primer sequences leads to a SyntaxError
igdiscover init --db=testdata/database --reads=testdata/reads.1.fastq.gz testrun/primers
pushd testrun/primers
igdiscover config \
	--set forward_primers "['CGTGA']" \
	--set reverse_primers "['TTCAC']"
igdiscover run -n stats/reads.json
popd

# Test using FLASH and parsing its log output
igdiscover init --db=testdata/database --reads=testdata/reads.1.fastq.gz testrun/flash
pushd testrun/flash
igdiscover config --set merge_program flash
igdiscover run stats/reads.json
popd


igdiscover init --db=testdata/database --reads=testdata/reads.1.fastq.gz testrun/paired
pushd testrun/paired
igdiscover config --set barcode_length_3prime 21

igdiscover run nofinal
if [[ -d final/ ]]; then
	echo "ERROR: nofinal failed"
	exit 1
fi

# run final iteration
igdiscover run

igdiscover run iteration-01/exact.tab
popd

# Use the merged file from above as input again
igdiscover init --db=testdata/database --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefastq && igdiscover run stats/reads.json )

# Test FASTA input
cutadapt --quiet -o testrun/reads.fasta testrun/paired/reads/2-merged.fastq.gz
igdiscover init --db=testdata/database --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefasta && igdiscover run stats/reads.json )
