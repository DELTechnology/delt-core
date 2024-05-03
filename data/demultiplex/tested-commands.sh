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