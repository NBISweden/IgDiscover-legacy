#!/bin/bash
set -eo pipefail
unset DISPLAY
export IGDATA=$PWD/igblast
source venv/bin/activate
rm -rf testrun
igdiscover init --library-name testing --db igdiscover-test-0.2/db --reads igdiscover-test-0.2/reads.1.fastq.gz testrun
cp -p igdiscover-test-0.2/igdiscover.yaml testrun/
cd testrun
igdiscover run -p
