#!/bin/bash
set -eo pipefail
unset DISPLAY
export IGDATA=$PWD/igblast
source venv/bin/activate
rm -rf testrun
igdiscover init --library-name=testing --db=testdata/db --reads testdata/reads.1.fastq.gz testrun
cp -p testdata/igdiscover.yaml testrun/
cd testrun
igdiscover run -p
