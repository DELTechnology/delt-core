#!/bin/bash
# make sure you installed pigz with `brew install pigz` to enable parallel processing

mkdir "eval"
ln -sf "../../../../../../../../polybox - Adriano Martinelli (adriano.martinelli@pharma.ethz.ch)@polybox.ethz.ch/decl-data/raw-files/o306001_1-sample1_704_504_OST1_S55_R1_001.fastq.gz" "eval/out.fastq.gz"

mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/0.S.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/0.S.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/0.S.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/1.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/1.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/1.C.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/2.B.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/2.B.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/2.B.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/3.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/3.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/3.C.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/4.B.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/4.B.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/4.B.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/5.C.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/5.C.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/5.C.cutadapt.log'


mv eval/out.fastq.gz eval/input.fastq.gz

cutadapt "eval/input.fastq.gz" \
-o "eval/out.fastq.gz" \
-e 0 \
-g "^file:cutadapt_input_files/6.S.fastq" \
--rename '{id} {comment}?{adapter_name}' \
--discard-untrimmed \
--json='/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/6.S.cutadapt.json' \
--cores=11 2>&1 | tee '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0/eval/6.S.cutadapt.log'

zgrep "@" "eval/out.fastq.gz" | gzip -c > "reads_with_adapters.gz"

python /Users/adrianomartinelli/projects/delt/delt-core/data/demultiplex/postprocess_output.py