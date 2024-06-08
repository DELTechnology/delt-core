import gzip
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import json

# Path to the gzip file
gzip_file_path = Path('eval/reads_with_adapters.gz').resolve()
assert gzip_file_path.exists()

path_to_dir = gzip_file_path.parent
number_of_reads = json.load(open(sorted((path_to_dir).glob('*.cutadapt.json'))[-1]))['read_counts']['output']

# Open the gzip file
from collections import defaultdict, Counter
selections = defaultdict(lambda: defaultdict(int))
with gzip.open(gzip_file_path, 'rt') as f:
    # Read the file line by line
    for line in tqdm(f, total=number_of_reads, ncols=100):
        # Split the line by '-'
        _, *adapters = line.strip().split('?')
        selection_id = '_'.join([i.split('.')[-1] for i in filter(lambda x: 'S' in x, adapters)])
        barcodes = tuple(i.split('.')[-1] for i in filter(lambda x: 'B' in x, adapters))

        selections[selection_id][barcodes] += 1

for selection_id, counts in tqdm(selections.items(), ncols=100):
    # counts = [(*i, j) for i,j in zip(counts.keys(), counts.values())]
    counts = [(j, *i) for i,j in zip(counts.keys(), counts.values())]
    # df = pd.DataFrame.from_records(counts, columns=['barcode1', 'barcode2', 'count'])
    df = pd.DataFrame.from_records(counts, columns=['count', 'barcode1', 'barcode2'])
    df = df.astype(int)
    df.sort_values(['barcode1', 'barcode2'], inplace=True)

    path_counts = path_to_dir / f'{selection_id}.tsv'
    df.to_csv(path_counts, index=False, sep='\t')
