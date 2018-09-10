#!/bin/bash
# Run this within an activated igdiscover environment
set -euo pipefail
unset DISPLAY

pytest

rm -rf testrun
mkdir testrun
[[ -L testdata ]] || ln -s igdiscover-testdata testdata

igdiscover init --db=testdata/database --reads=testdata/reads.1.fastq.gz testrun/paired

cd testrun/paired
igdiscover config --set barcode_length_3prime 21
igdiscover run
igdiscover run iteration-01/exact.tab
cd ../..

# Use the merged file from above as input again
igdiscover init --db=testdata/database --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefastq && igdiscover run stats/reads.json )

# Test FASTA input
sqt fastxmod -w 0 --fasta testrun/paired/reads/2-merged.fastq.gz > testrun/reads.fasta
igdiscover init --db=testdata/database --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefasta && igdiscover run stats/reads.json )
