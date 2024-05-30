from pathlib import Path

from pydantic import BaseModel, computed_field

path_to_struct = Path(
    '/Users/adrianomartinelli/polybox - Adriano Martinelli (adriano.martinelli@pharma.ethz.ch)@polybox.ethz.ch/decl-data/raw-files-downsampled/structure.txt')
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
import stat
import multiprocessing

path_demultiplex_exec = path_save / 'demultiplex.sh'

# note: delete existing demultiplex.sh
path_demultiplex_exec.unlink(missing_ok=True)

path_tmp_input_fastq = path_output_dir / 'input.fastq.gz'
path_output_fastq = path_output_dir / 'out.fastq.gz'
with open(path_save / 'demultiplex.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# make sure you installed pigz with `brew install pigz` to enable parallel processing\n\n')
    f.write(f'mkdir "{path_output_dir.relative_to(path_to_dir)}"\n')
    # we symlink the fastq file we want to demultiplex
    f.write(f'ln -sf "{os.path.relpath(path_input_fastq, start=path_output_dir)}" "{path_output_fastq.relative_to(path_to_dir)}"\n')

n_cores = multiprocessing.cpu_count()

error_rate = 1
rename_command = '{id} {comment}?{adapter_name}'

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

    report_file_name = f'{region.region_id}.cutadapt.json'
    stdout_file_name = f'{region.region_id}.cutadapt.log'
    info_file_name = f'{region.region_id}.cutadapt.info.gz'

    with open(path_demultiplex_exec, 'a') as f:
        cmd = f"""
            mv {path_output_fastq} {path_input_fastq}
            
            cutadapt "{path_input_fastq}" \\
            -o "{path_output_fastq}" \\
            -e {error_rate} \\
            -g "^file:{path_adapters}" \\
            --rename '{rename_command}' \\
            --discard-untrimmed \\
            --json={report_file_name} \\
            --info-file={info_file_name} \\
            --cores={n_cores} 2>&1 | tee "{stdout_file_name}"
            
            """
        cmd = textwrap.dedent(cmd)
        f.write(cmd)

with open(path_demultiplex_exec, 'a') as f:
    f.write(f'zgrep "@" "{path_output_fastq}" | gzip -c > "reads_with_adapters.gz"')

os.chmod(path_demultiplex_exec, os.stat(path_demultiplex_exec).st_mode | stat.S_IEXEC)
