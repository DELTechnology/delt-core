import pandas as pd
from pathlib import Path

sel_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/experiments/test-1/selections')

neg_ctrls = ['AG24_1', 'AG24_2', 'AG24_3']
naive = ['AG24_19', 'AG24_20', 'AG24_21']
protein = ['AG24_10', 'AG24_11', 'AG24_12']

data = []
for name in neg_ctrls:
    df = pd.read_csv(sel_dir / name / 'counts.txt', sep='\t')
    df['name'] = name
    df['group'] = 'neg_ctrl'
    data.append(df)

for name in naive:
    df = pd.read_csv(sel_dir / name / 'counts.txt', sep='\t')
    df['name'] = name
    df['group'] = 'naive'
    data.append(df)

for name in protein:
    df = pd.read_csv(sel_dir / name / 'counts.txt', sep='\t')
    df['name'] = name
    df['group'] = 'protein'
    data.append(df)

data = pd.concat(data)

data.groupby('name')['count'].sum()


