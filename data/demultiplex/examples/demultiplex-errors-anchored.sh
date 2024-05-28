# sequential demultiplexing with exact matches
# we allow one error, e > 1/3
cutadapt -e 1 \
-g ^AAA \
-o out.fastq input.fastq \
--rename '{id}-{comment}-{adapter_name}' \
--info-file=info.tsv

cutadapt -e 2 \
-g "^AAA" \
-o out.fastq input.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.tsv

# NOTE: if we allow errors and we have indels at the end of the sequences
# cutadapt favours deletions over insertions, i.e.
# for an adapter ^AAA we see

# A	TAABTCAAABTCAAA
# AA	TABTCAAABTCAAA
