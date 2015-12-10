#!/bin/bash
set -eo pipefail
unset DISPLAY

export IGDATA=$HOME/scm/antibodies/igblast
source venv/bin/activate
rm -rf tmp/testrun
igypipe init --library-name testing --db testdata/db --reads testdata/test.1.fastq.gz tmp/testrun
cp -p testdata/igypipe.yaml tmp/testrun/
cd tmp/testrun
snakemake -pj
