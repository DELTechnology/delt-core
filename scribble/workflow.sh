#!/usr/bin/env bash
# excel_path=Path('/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx)
# config_path=Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml')

delt-cli demultiplex init --excel_path=/Users/adrianomartinelli/projects/delt/delt-core/templates/NF.xlsx
delt-cli demultiplex run --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml

delt-cli analyse add \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test-1 \
--selections='["AG24_1", "AG24_2", "AG24_3", "AG24_10", "AG24_11", "AG24_12", "AG24_19", "AG24_20", "AG24_21"]'

delt-cli analyse add \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test-2 \
--selections='[AG24_28, "AG24_29", "AG24_30", "AG24_37", "AG24_38", "AG24_39", "AG24_46", "AG24_47", "AG24_48"]'

delt-cli analyse enrichment \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test-1 \
--method=counts

delt-cli analyse enrichment \
--config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
--name=test-1 \
--method=edgeR

delt-cli library enumerate --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml
delt-cli library properties --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml
delt-cli library represent --method=morgan --config_path=/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/config.yaml \
