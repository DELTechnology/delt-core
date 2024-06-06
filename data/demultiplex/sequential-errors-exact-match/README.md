# Potential Optimizations
## Considerations for Large File Sets
Memory Usage: Running awk on very large files or many files at once can be resource-intensive.
Monitor system resources if this becomes a concern.
Parallel Processing: If you're working with an extremely large number of files or very large files, consider parallel 
processing techniques to speed up processing. Tools like parallel can run multiple awk instances concurrently:
```bash
parallel "awk -F'\t' 'BEGIN {{OFS=\"\t\"}} $2 != \"-1\" {{print \$1, \$8}}' {} > {.}_filtered.tsv" ::: info.*.tsv
```
This uses GNU parallel to run an awk command on each file simultaneously, outputting to individual files for each input.