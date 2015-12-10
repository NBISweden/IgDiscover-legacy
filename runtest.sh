#!/bin/bash
set -eo pipefail
unset DISPLAY
export IGDATA=$PWD/igblast
source venv/bin/activate
rm -rf testrun
igypipe init --library-name testing --db testdata/db --reads testdata/reads.1.fastq.gz testrun
cp -p testdata/igypipe.yaml testrun/
cd testrun
snakemake -pj
