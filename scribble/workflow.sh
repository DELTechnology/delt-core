#!/usr/bin/env bash
# excel_path=Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx)
# config_path=Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')

delt-cli demultiplex init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx
delt-cli demultiplex run --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml

delt-cli analyse add \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test \
--selections='["AG24_1", "AG24_2", "AG24_3", "AG24_10", "AG24_11", "AG24_12", "AG24_19", "AG24_20", "AG24_21"]'

delt-cli analyse enrichment \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test \
--method=counts

delt-cli analyse enrichment \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test \
--method=edgeR

delt-cli properties run --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml

delt-cli represent run --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--method=morgan