# clusterplot test data

The `clusterplot.tab.gz` file was extracted from the SRR5408020 dataset with this command:

    (
      zcat iteration-01/filtered.tab.gz | head -n 1
      zcat iteration-01/filtered.tab.gz | awk -vFS='\t' '$19==0' | grep '\bIGHV1-18\*01\b' | head -n 10
    ) | cut -f2,19,38 | gzip > clusterplot.tab

# IgBLAST parsing

igblast.fasta is a test file for IgBLAST parsing with some corner cases (V incomplete at 5'
end, no D assigned etc.). The database for this is in database/.


