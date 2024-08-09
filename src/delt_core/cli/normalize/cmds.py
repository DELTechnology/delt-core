from pathlib import Path
from importlib import resources
import subprocess

from ... import normalize as n
from delt_core.demultiplex.validation import Config


def run(
        config_file: Path,
        target = list,
        control = list,
) -> None:
    config = Config.from_yaml(config_file).model_dump()
    
    root = config['Root']
    selection_file = config['Selection']['SelectionFile']
    data_dir = 'evaluations'
    output_dir = 'normalization'
    target_ids = ','.join(target.split())
    control_ids = ','.join(control.split())

    args = [f'{root} {selection_file} {data_dir} {output_dir} {target_ids} {control_ids}']
    Path(output_dir).mkdir(exist_ok=True)

    script = 'normalization.R'
    dir = 'delt_core.normalize'
    with resources.path(dir, script) as path:
        cmd = f'Rscript {str(path)}'
        subprocess.run(cmd.split() + args, check=True)

