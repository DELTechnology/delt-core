from pathlib import Path

from pydantic import BaseModel, computed_field

path_to_struct = Path(
    '/Users/adrianomartinelli/polybox - Adriano Martinelli (adriano.martinelli@pharma.ethz.ch)@polybox.ethz.ch/decl-data/raw-files-downsampled/structureNF2.txt')
path_to_dir = path_to_struct.parent
path_save = path_to_dir / 'cutadapt_input_files'
path_save.mkdir(parents=True, exist_ok=True)

with open(path_to_struct, 'r') as f:
    struct = f.read().split('\n')

path_input_fastq = path_to_dir / Path(struct.pop(0))
path_input_fastq = path_input_fastq.resolve()

path_output_dir = path_to_dir / struct.pop(0)
path_output_dir = path_output_dir.resolve()

class Region(BaseModel):
    start: int
    end: int
    type: str
    codons: list[str]
    position_in_construct: int = None

    @computed_field
    @property
    def region_id(self) -> str:
        return f'{self.position_in_construct}.{self.type}'


regions = []
for i in struct:
    if len(i) == 0:
        continue
    idx_start, idx_end, _type, relative_path_to_codons = i.split('\t')

    with open(path_to_dir / relative_path_to_codons, 'r') as f:
        codons = f.read().split('\n')
        codons = filter(len, codons)

    region = Region(start=int(idx_start) - 1, end=int(idx_end) - 1, type=_type, codons=codons)
    regions.append(region)

regions = sorted(regions, key=lambda x: x.start)

fastq = ''
for i, region in enumerate(regions):
    region.position_in_construct = i
    fastq = [f'>{region.region_id}.{index}\n{codon}'
             for index, codon in enumerate(region.codons)]
    fastq = '\n'.join(fastq)

    with open(path_save / f'{region.region_id}.fastq', 'w') as f:
        f.write(fastq)

# %%
import textwrap
import os

# note: delete existing demultiplex.sh
(path_save / 'demultiplex.sh').unlink(missing_ok=True)

path_tmp_input_fastq = path_output_dir / 'input.fastq.gz'
with open(path_save / 'demultiplex.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
    f.write(f'mkdir "{path_output_dir.relative_to(path_to_dir)}"\n')
    # we symlink the fastq file we want to demultiplex
    f.write(f'ln -sf "{os.path.relpath(path_input_fastq, start=path_output_dir)}" "{path_tmp_input_fastq}"\n')

import multiprocessing
n_cores = multiprocessing.cpu_count()

error_rate = 0
rename_command = '{id}-{comment}-{adapter_name}'

use_relative_paths = True
for region in regions:
    path_output_fastq = path_output_dir / 'out.fastq.gz'
    # path_output_info = path_output_dir / 'info.tsv'
    path_adapters = path_save / f'{region.region_id}.fastq'
    # note: from now on we use the output of the previous step as input
    path_input_fastq = path_output_dir / 'input.fastq.gz'

    if use_relative_paths:
        path_adapters = path_adapters.absolute().relative_to(path_to_dir)
        path_output_fastq = path_output_fastq.absolute().relative_to(path_to_dir)
        # path_output_info = path_output_info.absolute().relative_to(path_to_dir)
        path_input_fastq = path_input_fastq.absolute().relative_to(path_to_dir)

    with open(path_save / 'demultiplex.sh', 'a') as f:
        cmd = f"""
            cutadapt "{path_input_fastq}" \\
            -o "{path_output_fastq}" \\
            -e {error_rate} \\
            -g "^file:{path_adapters}" \\
            --rename '{rename_command}' \\
            --discard-untrimmed \\
            --cores={n_cores}
            
            mv {path_output_fastq} {path_input_fastq}
            """
        cmd = textwrap.dedent(cmd)
        f.write(cmd)


# cutadapt -e 0 \
# -g "^file:cutadapt_input_files/0.S.fastq" \
# -o "Evaluation_2301_704504_NF2-AG_yOST/out.fastq.gz" "o306001_1-sample1_704_504_OST1_S55_R1_001.fastq.gz" \
# --rename '{id}-{comment}-{adapter_name}' \
# --discard-untrimmed \
# --cores=5
#
# cutadapt -e 0 \
# -g "^file:cutadapt_input_files/0.S.fastq" \
# -o "Evaluation_2301_704504_NF2-AG_yOST/out.fastq.gz" "../raw-files/o306001_1-sample1_704_504_OST1_S55_R1_001.fastq.gz" \
# --rename '{id}-{comment}-{adapter_name}' \
# --discard-untrimmed \
# --cores=5