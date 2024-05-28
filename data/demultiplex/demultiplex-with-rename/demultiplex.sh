# sequential demultiplexing with exact matches
cutadapt -e 0 \
-g "file:../p1.fasta;min_overlap=6" \
-o out.fastq ../simlation.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.p1.tsv \
--discard-untrimmed

cutadapt -e 0 \
-g "file:../c1.fasta;min_overlap=24" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.c1.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

cutadapt -e 0 \
-g "file:../b1.fasta;min_overlap=6" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.b1.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

cutadapt -e 0 \
-g "file:../c2.fasta;min_overlap=25" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.c2.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

cutadapt -e 0 \
-g "file:../b2.fasta;min_overlap=7" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.b2.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

cutadapt -e 0 \
-g "file:../c3.fasta;min_overlap=12" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.c3.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

cutadapt -e 0 \
-g "file:../p2.fasta;min_overlap=9" \
-o tmp.fastq out.fastq \
--rename '{id} {comment} {adapter_name}' \
--info-file=info.p2.tsv \
--discard-untrimmed && mv tmp.fastq out.fastq

# get all read names
grep @ out.fastq > 'valid_reads.txt'

# remove unmatched reads from info files
# for i in info.*.tsv;
#do
#  grep -v '\t-1\t' $i > temp_file && mv temp_file $i
#done

# extract adapter names
#awk -F'\t' 'BEGIN {OFS="\t"} $2 != "-1" {print $1, $8}' info.p1.tsv
#awk -F'\t' 'BEGIN {OFS="\t"} $2 != "-1" {print $1, $8}' info.*.tsv > info.tsv

