#!/usr/bin/env python3
from cutadapt import seqio
import sys
from collections import defaultdict

with seqio.open(sys.argv[1]) as sr:
	reads = list(sr)

d = defaultdict(list)
for read in reads:
	d[read.sequence].append(read.name)

for seq, names in d.items():
	print('>{}\n{}'.format(','.join(names), seq))

