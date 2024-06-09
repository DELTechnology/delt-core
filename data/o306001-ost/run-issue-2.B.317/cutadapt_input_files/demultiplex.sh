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
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/0.S.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/0.S.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/1.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/1.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/1.C.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.2.B.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/2.B.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/2.B.cutadapt.json' \
--info-file='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/2.B.cutadapt.info.gz' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/2.B.cutadapt.log'

ln -sf out.2.B.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.3.C.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/3.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--untrimmed-output "eval/3.C.untrimmed.output.gz" \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/3.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/3.C.cutadapt.log'


ln -sf out.3.C.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.4.B.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/4.B.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--untrimmed-output "eval/4.B.untrimmed.output.gz" \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/4.B.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/4.B.cutadapt.log'


ln -sf out.4.B.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.5.C.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/5.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--untrimmed-output "eval/5.C.untrimmed.output.gz" \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/5.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/5.C.cutadapt.log'


ln -sf out.5.C.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.6.S.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/6.S.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--untrimmed-output "eval/6.S.untrimmed.output.gz" \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/6.S.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-issue-2.B.317/eval/6.S.cutadapt.log'

#zgrep "@" "eval/out.fastq.gz" | gzip -c > "eval/reads_with_adapters.gz"
#rm "eval/out.fastq.gz" "eval/input.fastq.gz"