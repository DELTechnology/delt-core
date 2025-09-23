#!/usr/bin/env bash
# excel_path=Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx)
# config_path=Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')

delt-cli demultiplex init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx
delt-cli demultiplex run --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml

delt-cli analyze enrichment --method=counts --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml
delt-cli analyze enrichment --method=edgeR
#delt-cli analyze enrichment --method=DESeq2