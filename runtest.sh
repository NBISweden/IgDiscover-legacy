#!/bin/bash
# Run this within an activated igdiscover environment

set -eo pipefail
unset DISPLAY

rm -rf testrun
mkdir testrun
igdiscover init --db=testdata/db --reads=testdata/reads.1.fastq.gz testrun/paired

cd testrun/paired
igdiscover config --set iterations 1 \
   --set race_g true \
   --set stranded true \
   --set barcode_length_5prime 12 \
   --set cdr3_location '[-80, -60]' \
   --set reverse_primers '["GCAGGCCTTTTTGGCCNNNNNGGGGCATTCTCACAGGAGACGAGGGGGAAAAG"]' \
   --set forward_primers '["CGTGAGCTGAGTACGACTCACTATAGCTTCAC"]' \
   --set germline_filter.unique_js 1
time -p igdiscover run -p
cd ../..

# Use the merged file from above as input again
igdiscover init --db=testdata/db --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefastq && exec time -p igdiscover run -p stats/reads.json )

# Test FASTA input
sqt fastxmod -w 0 --fasta testrun/paired/reads/2-merged.fastq.gz > testrun/reads.fasta
igdiscover init --db=testdata/db --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefasta && exec time -p igdiscover run -p stats/reads.json )
