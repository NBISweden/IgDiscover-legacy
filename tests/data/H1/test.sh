#!/bin/bash
set -euo pipefail

igdiscover germlinefilter --whitelist=V.fasta --max-differences=0 --unique-CDR3=5 --cluster-size=100 --unique-J=3 --cross-mapping-ratio=0.02 --allele-ratio=0.1 candidates.tab.gz > new_V_germline.tab
diff -U 0 <(cut -f1 expected.tab) <(cut -f1 new_V_germline.tab) | grep '^[+-]' | sed 1,2d
