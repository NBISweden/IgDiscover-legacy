#!/bin/bash
set -eo pipefail
unset DISPLAY
export IGDATA=$PWD/igblast
source venv/bin/activate

rm -rf testrun
mkdir testrun
igdiscover init --library-name=testing --db=testdata/db --reads=testdata/reads.1.fastq.gz testrun/paired
cp -p testdata/igdiscover.yaml testrun/paired/
cd testrun/paired
time igdiscover run -p
cd ../..

# Use the merged file from above as input again
igdiscover init --library-name=testing --db=testdata/db --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testdata/igdiscover.yaml testrun/singlefastq/
cd testrun/singlefastq
time igdiscover run -p stats/reads.json
cd ../..

sqt fastxmod -w 0 --fasta testrun/paired/reads/2-merged.fastq.gz > testrun/reads.fasta
igdiscover init --library-name=testing --db=testdata/db --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testdata/igdiscover.yaml testrun/singlefasta/
cd testrun/singlefasta
time igdiscover run -p stats/reads.json
