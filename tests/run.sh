#!/bin/bash
# Run this within an activated igdiscover environment
set -euo pipefail
unset DISPLAY

nosetests --with-doctest -P tests/ igdiscover/

rm -rf testrun
mkdir testrun
[[ -L testdata ]] || ln -s igdiscover-testdata-* testdata

igdiscover init --db=testdata/db --reads=testdata/reads.1.fastq.gz testrun/paired

cd testrun/paired
igdiscover config --set iterations 2 \
   --set race_g true \
   --set stranded true \
   --set barcode_length_5prime 12 \
   --set cdr3_location '[-80, -60]' \
   --set reverse_primers '["GCAGGCCTTTTTGGCCNNNNNGGGGCATTCTCACAGGAGACGAGGGGGAAAAG"]' \
   --set forward_primers '["CGTGAGCTGAGTACGACTCACTATAGCTTCAC"]' \
   --set germline_filter.unique_js 1 \
   --set j_discovery.propagate false
igdiscover run
cd ../..

# Use the merged file from above as input again
igdiscover init --db=testdata/db --single-reads=testrun/paired/reads/2-merged.fastq.gz testrun/singlefastq
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefastq && igdiscover run stats/reads.json )

# Test FASTA input
sqt fastxmod -w 0 --fasta testrun/paired/reads/2-merged.fastq.gz > testrun/reads.fasta
igdiscover init --db=testdata/db --single-reads=testrun/reads.fasta testrun/singlefasta
cp -p testrun/paired/igdiscover.yaml testrun/singlefastq/
( cd testrun/singlefasta && igdiscover run stats/reads.json )
