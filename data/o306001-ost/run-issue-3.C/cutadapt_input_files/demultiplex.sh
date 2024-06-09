#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "eval"
ln -sf "../../../../../../../../polybox - Adriano Martinelli (adriano.martinelli@pharma.ethz.ch)@polybox.ethz.ch/decl-data/raw-files/o306001_1-sample1_704_504_OST1_S55_R1_001.fastq.gz" "eval/out.fastq.gz"

mv eval/out.fastq.gz eval/input.fastq.gz

gzcat eval/input.fastq.gz | head -4000000 | cutadapt - \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/0.S.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='eval/0.S.cutadapt.json' \
--cores=11 2>&1 | tee 'eval/0.S.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/1.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='eval/1.C.cutadapt.json' \
--cores=11 2>&1 | tee 'eval/1.C.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/2.B.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='eval/2.B.cutadapt.json' \
--cores=11 2>&1 | tee 'eval/2.B.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/3.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--untrimmed-output eval/untrimmed.gz \
--info-file='eval/3.C.info.gz' \
--json='eval/3.C.cutadapt.json' \
--cores=11 2>&1 | tee 'eval/3.C.cutadapt.log'


#gzcat input.fastq.gz | head -1000 | \
#cutadapt -a CGAGTCCCATGGCGCCGGATCGACG \
#--untrimmed-output eval/untrimmed.gz -o eval/trimmed-output.gz --info-file eval/info.tsv -
