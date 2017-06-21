#!/bin/bash
# Run this within an activated igdiscover environment

set -eo pipefail
unset DISPLAY

rm -rf testrun
mkdir testrun
igdiscover init --db=testdata/db --reads=testdata/reads.1.fastq.gz testrun/paired
cp -p testdata/igdiscover.yaml testrun/paired/
( cd testrun/paired && exec time -p igdiscover run -p )

# Use the merged file from above as input again
igdiscover init --db=testdata/db --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testdata/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefastq && exec time -p igdiscover run -p stats/reads.json )

sqt fastxmod -w 0 --fasta testrun/paired/reads/2-merged.fastq.gz > testrun/reads.fasta
igdiscover init --db=testdata/db --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testdata/igdiscover.yaml testrun/singlefasta/
( cd testrun/singlefasta && exec time -p igdiscover run -p stats/reads.json )
