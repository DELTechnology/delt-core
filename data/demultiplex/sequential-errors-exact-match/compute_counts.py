from collections import defaultdict
from pathlib import Path

reads = defaultdict(list)

path_valid_read_names = '/Users/adrianomartinelli/projects/delt/delt-core/data/demultiplex/sequential-errors-exact-match/valid_reads.txt'
with open(path_valid_read_names, 'r') as f:
    valid_read_names = {i.strip('@').strip('\n') for i in f.readlines()}

fname = 'info.tsv'
fname = '/Users/adrianomartinelli/projects/delt/delt-core/data/demultiplex/sequential-errors-exact-match/info.tsv'
with open(fname, 'r') as f:
    for line in f:
        read_name, adapter_id = line.strip('\n').split('\t')
        if read_name not in valid_read_names:
            continue
        if adapter_id.startswith('c'):
            continue
        pos, idx = adapter_id.split('.')
        reads[read_name].append((pos, int(idx)))

counts = defaultdict(lambda: defaultdict(int))
for read_name, codes in reads.items():
    bbs = tuple([i[1] for i in sorted(filter(lambda x: x[0].startswith('bb'), codes), key=lambda x: x[1])])
    ps = tuple([i[1] for i in sorted(filter(lambda x: x[0].startswith('p'), codes), key=lambda x: x[1])])
    if bbs:
        counts[ps][bbs] += 1

n_bbs = len(bbs)

parent_dir = Path(fname).parent
for sel_id, sel_counts in counts.items():
    sel_id = '_'.join([str(i) for i in sel_id])
    fpath = parent_dir / f'counts_{sel_id}.tsv'

    with open(fpath, 'w') as f:
        headers = '\t'.join([f'bb{i+1}' for i in range(n_bbs)] + ['counts\n'])
        f.write(headers)

        for bbs, count in sel_counts.items():
            bbs = '\t'.join([str(i) for i in bbs])
            f.write(f'{bbs}\t{count}\n')
