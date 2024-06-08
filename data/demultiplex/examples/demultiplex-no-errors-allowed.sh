# sequential demultiplexing with exact matches
cutadapt -e 0 \
-g "AAA;min_overlap=3" \
-o out.fastq input.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.tsv
