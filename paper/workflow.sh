#!/usr/bin/env bash
# excel_path=Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx)
# config_path=Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config_with_steps.yaml')

delt-hit init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/NF.xlsx
delt-hit demultiplex prepare --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml
paper/experiment-1/demultiplex/cutadapt_input_files/demultiplex.sh

delt-hit demultiplex report --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml
delt-hit demultiplex qc --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml

delt-hit analyse enrichment \
--config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml \
--name=analysis-1 \
--method=counts

delt-hit analyse enrichment \
--config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml \
--name=analysis-1 \
--method=edgeR

delt-hit analyse enrichment \
--config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml \
--name=analysis-2 \
--method=counts

delt-hit analyse enrichment \
--config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml \
--name=analysis-2 \
--method=edgeR

delt-hit library enumerate --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml
delt-hit library properties --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml
delt-hit library represent --method=morgan --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml
delt-hit library represent --method=bert --config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml

delt-hit dashboard \
--config_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/config.yaml \
--counts_path=/Users/adrianomartinelli/projects/delt/delt-core/paper/experiment-1/selections/AG24_10/counts.txt
