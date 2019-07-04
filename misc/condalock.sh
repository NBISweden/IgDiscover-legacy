#!/bin/bash
set -xeuo pipefail
conda env create -n igdlock -f environment.yml
conda env export -n igdlock | grep -Ev '^(name|prefix):' > environment.lock.yml
conda env remove -n igdlock
