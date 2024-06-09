# Run Analysis
Symlink the fastq file into the different `run-*` folders

```bash
for d in ls -d run-*; do
    ln -s <PATH_TO_FASTQ> $d/input.fastq.gz
done
```
Run `demultiplex.sh` in each `run-*` folder
```bash
for d in ls -d run-*; do
    cd $d
    bash cutadapt_input_files/demultiplex.sh
    cd ..
done
``` 
