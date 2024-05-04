# most basic command
cutadapt -g GCCTCG -o out.fastq simlation.fastq --info-file=info.tsv

# use fasta file to specify adapters
cutadapt -g file:c1.fasta -o out.fastq simlation.fastq --info-file=info.tsv

# trim multiple adapter sequences
cutadapt -g ACACGGAGCTT -g CTGTGTGCTGGC \
-o out.fastq simlation.fastq \
--info-file=info.tsv

# trim multiple adapter sequences
cutadapt -g GGAGCTT -a CGAGTCCC -n 2 \
-o out.fastq simlation.fastq \
--info-file=info.tsv

cutadapt \
-g GGAGCTTCTGAATTCTGTGTGCTG \
-g CGAGTCCCATGGCGCCGGATCGACG \
-g GCGTCAGGCAGC \
-n 3 \
-o out.fastq simlation.fastq \
--info-file=info.tsv

cutadapt \
-g file:p1.fasta \
-o out1.fastq simlation.fastq \
--info-file=info.tsv

cutadapt \
-g file:c1.fasta \
-o out2.fastq out1.fastq \
--info-file=info.tsv

cutadapt \
-g file:b1.fasta \
-o out3.fastq out2.fastq \
--info-file=info.tsv

cutadapt \
-g file:c2.fasta \
-o out4.fastq out3.fastq \
--info-file=info.tsv

cutadapt \
-g file:b2.fasta \
-o out5.fastq out4.fastq \
--info-file=info.tsv

cutadapt \
-g file:c3.fasta \
-o out6.fastq out5.fastq \
--info-file=info.tsv

cutadapt \
-g file:p2.fasta \
-o out7.fastq out6.fastq \
--info-file=info.tsv

# only match exact matches
cutadapt -e 0 -O 10 \
-g AAAAAAAAAA \
-o out.short.fastq short_reads.fastq \
--info-file=info.tsv

# allow 1 SNPs
cutadapt -e 0.1 -O 10 \
-g AAAAAAAAAA \
-o out.short.fastq short_reads.fastq \
--info-file=info.tsv

# allow 2 SNPs
cutadapt -e 0.2 -O 10 \
-g AAAAAAAAAA \
-o out.short.fastq short_reads.fastq \
--info-file=info.tsv

# allow del only
# however, this should allow for SNPs at the ends of the reads as well
# but I observe that it does only for SNPs at the end not the start of the read
cutadapt -e 0 -O 9 \
-g AAAAAAAAAA \
-o /dev/null short_reads.fastq \
--info-file=info.tsv

cutadapt -e 0 -O 9 \
-g AAAAAAAAAA \
-o /dev/null input.fastq \
--info-file=info.tsv

# demultiplex in different files
cutadapt -e 0 -O 6 \
-g file:p1.fasta \
-a file:p2.fasta \
-o out-{name}.fastq simlation.fastq \
--info-file=info.tsv

# sequential demultiplexing with exact matches
cutadapt -e 0 \
-g "file:p1.fasta;min_overlap=6" \
-o out1.fastq simlation.fastq \
--discard-untrimmed