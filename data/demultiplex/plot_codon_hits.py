from pathlib import Path
import json

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

path_exp_dir = Path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/delt/run-e0')
path_eval_dir = path_exp_dir / 'eval'
report_paths = path_eval_dir.glob('*.cutadapt.json')

cont = []
for report_path in report_paths:
    report = json.load(report_path.open('r'))
    stats = report['adapters_read1']
    for stat in stats:
        if stat['five_prime_end']['trimmed_lengths']:
            item = {}
            item['region_id'] = '_'.join(stat['name'].split('.')[:-1])
            item['index'] = int(stat['name'].split('.')[-1])
            item['expected'] = stat['five_prime_end']['trimmed_lengths'][0]['expect']
            item['counts'] = stat['five_prime_end']['trimmed_lengths'][0]['counts']
            item['counts'] = dict(zip(range(len(item['counts'])), item['counts']))
            cont.append(item)

df = pd.DataFrame(cont)
errors = df.counts.apply(pd.Series)
df = pd.concat([df, errors], axis=1).melt(id_vars=df.columns, var_name='number_of_errors', value_name='error_counts')

for grp_name, grp_dat in df.groupby('region_id'):
    # grp_dat = pd.concat([grp_dat, grp_dat.assign(number_of_errors=1)])
    pdat = grp_dat.pivot(index='index', columns='number_of_errors', values='error_counts')

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle(grp_name)

    pdat.plot(kind='bar', stacked=True, ax=axs[0])
    pdat.plot(kind='bar', stacked=True, ax=axs[1])

    axs[0].set_ylabel('counts')
    axs[1].set_yscale('log')

    expected = grp_dat.expected.iloc[0]
    axs[1].plot([0, grp_dat.index.max()], [expected, expected], 'k--')

    # g = sns.barplot(data=grp_dat,
    #                 x='index', y='error_counts',
    #                 hue='number_of_errors',
    #                 estimator='sum',
    #                 ax=axs[0])
    # g = sns.barplot(data=grp_dat,
    #                     x='index', y='error_counts',
    #                     hue='number_of_errors',
    #                 estimator='sum',
    #                     ax=axs[1])

    for ax in axs:
        ax.set_xticklabels([])

    fig.tight_layout()
    fig.savefig(path_eval_dir / f'hits_{grp_name}.pdf')
    fig.show()
