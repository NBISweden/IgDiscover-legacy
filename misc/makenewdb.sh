#!/bin/bash
set -euo pipefail
grep -v '^$' "Indian Rhesus germline IGHV-SundlingSTM.vMar19 2014.txt" | \
	dos2unix | \
	awk '
		NR%2==1 {
			name=substr($1, 2);
			sub(/^VH/, "IGHV", name);
			sub(/\./, "-", name);
			name = name "*01";
			printf("MacacaMulatta\tIGH\tVariable\t%s\t-1\t", name)
		}
		NR%2==0 {
			printf("%s\n", toupper($1))
		}' > segment.new.txt
(
	awk '!($1 == "MacacaMulatta" && $2 == "IGH" && $3 == "Variable")' macaca.imgt.txt
	cat segment.new.txt
) > segments.all.minor.txt
echo "new segments.all.minor.txt created"
./preparedb.py segments.all.minor.txt
echo "database/ created"
rm segment.new.txt
