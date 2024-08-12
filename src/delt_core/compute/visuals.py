from collections import defaultdict
import os
from pathlib import Path

from faerun import Faerun
from matplotlib.colors import ListedColormap
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.AllChem import GetRDKitFPGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
import scipy.stats as ss
import tmap as tm

from lipinski import lipinski_pass
from utils import read_txt


def read_counts(
        file: Path,
) -> dict:
    counts = defaultdict(int)
    lines = read_txt(file)[1:]
    for line in lines:
        count, *bbs = line.strip().split('\t')
        counts['-'.join(bbs)] = int(count)
    return counts


def compute_fingerprint(
        mol: Chem.rdchem.Mol,
        num_bits: int = 2048,
) -> ExplicitBitVect:
    fpgen = GetRDKitFPGenerator(fpSize=num_bits)
    return fpgen.GetFingerprint(mol)


def compute_tree(
        smiles: list,
        labels: list,
        counts: list = None,
) -> None:

    enc = tm.Minhash(128)
    lf = tm.LSHForest(2048, 128)

    fps = []
    tpsa = []
    logp = []
    mw = []
    h_acceptors = []
    h_donors = []
    ring_count = []
    is_lipinski = []
    has_coc = []
    has_sa = []
    has_tz = []
    
    substruct_coc = AllChem.MolFromSmiles('COC')
    substruct_sa = AllChem.MolFromSmiles('NS(=O)=O')
    substruct_tz = AllChem.MolFromSmiles('N1N=NN=C1')

    mols = [AllChem.MolFromSmiles(m) for m in smiles]
    total = len(mols)
    if not counts:
        counts = [0] * total

    for i, mol in enumerate(mols):
        if i % 1000 == 0 and i > 0:
            print(f'{round(100 * (i / total))}% done ...')

        fps.append(tm.VectorUchar(compute_fingerprint(mol)))
        tpsa.append(Descriptors.TPSA(mol))
        logp.append(Descriptors.MolLogP(mol))
        mw.append(Descriptors.MolWt(mol))
        h_acceptors.append(Descriptors.NumHAcceptors(mol))
        h_donors.append(Descriptors.NumHDonors(mol))
        ring_count.append(Descriptors.RingCount(mol))
        is_lipinski.append(lipinski_pass(mol))
        has_coc.append(mol.HasSubstructMatch(substruct_coc))
        has_sa.append(mol.HasSubstructMatch(substruct_sa))
        has_tz.append(mol.HasSubstructMatch(substruct_tz))

    counts_ranked = ss.rankdata(np.array(counts) / max(counts)) / len(counts)
    tpsa_ranked = ss.rankdata(np.array(tpsa) / max(tpsa)) / len(tpsa)
    logp_ranked = ss.rankdata(np.array(logp) / max(logp)) / len(logp)
    mw_ranked = ss.rankdata(np.array(mw) / max(mw)) / len(mw)
    h_acceptors_ranked = ss.rankdata(np.array(h_acceptors) / max(h_acceptors)) / len(h_acceptors)
    h_donors_ranked = ss.rankdata(np.array(h_donors) / max(h_donors)) / len(h_donors)
    ring_count_ranked = ss.rankdata(np.array(ring_count) / max(ring_count)) / len(ring_count)

    fps = enc.batch_from_binary_array(fps)
    lf.batch_add(fps)
    lf.index()
    cfg = tm.LayoutConfiguration()
    cfg.k = 100
    # cfg.sl_extra_scaling_steps = 1
    cfg.sl_repeats = 2
    cfg.mmm_repeats = 2
    cfg.node_size = 2
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)

    bin_cmap = ListedColormap(['#e74c3c', '#2ecc71'], name='bin_cmap')

    f = Faerun(
        clear_color='#222222',
        coords=False,
        view='front',
        impress='made with <a href="http://tmap.gdb.tools" target="_blank">tmap</a><br />and <a href="https://github.com/reymond-group/faerun-python" target="_blank">faerun</a>',
    )

    f.add_scatter(
        'Property',
        {
            'x': x,
            'y': y,
            'c': [
                counts_ranked,
                is_lipinski,
                has_coc,
                has_sa,
                has_tz,
                tpsa_ranked,
                logp_ranked,
                mw_ranked,
                h_acceptors_ranked,
                h_donors_ranked,
                ring_count_ranked,
            ],
            'labels': labels,
        },
        shader='smoothCircle',
        colormap=[
            'viridis',
            bin_cmap,
            bin_cmap,
            bin_cmap,
            bin_cmap,
            'viridis',
            'viridis',
            'viridis',
            'viridis',
            'viridis',
            'viridis',
        ],
        point_scale=2.5,
        categorical=[False, True, True, True, True, False, False, False, False, False],
        has_legend=True,
        legend_labels=[
            None,
            [(0, 'No'), (1, 'Yes')],
            [(0, 'No'), (1, 'Yes')],
            [(0, 'No'), (1, 'Yes')],
            [(0, 'No'), (1, 'Yes')],
        ],
        selected_labels=['SMILES', 'ID'],
        series_title=[
            'Counts',
            'Lipinski',
            'Ethers',
            'Sulfonamides',
            'Tetrazoles',
            'TPSA',
            'logP',
            'Mol Weight',
            'H Acceptors',
            'H Donors',
            'Ring Count',
        ],
        max_legend_label=[
            str(round(max(counts))),
            None,
            None,
            None,
            None,
            str(round(max(tpsa))),
            str(round(max(logp))),
            str(round(max(mw))),
            str(round(max(h_acceptors))),
            str(round(max(h_donors))),
            str(round(max(ring_count))),
        ],
        min_legend_label=[
            str(round(min(counts))),
            None,
            None,
            None,
            None,
            str(round(min(tpsa))),
            str(round(min(logp))),
            str(round(min(mw))),
            str(round(min(h_acceptors))),
            str(round(min(h_donors))),
            str(round(min(ring_count))),
        ],
        title_index=2,
        legend_title='Property',
    )

    f.add_tree('deltree', {'from': s, 'to': t}, point_helper='Property')
    f.plot('tree')


def run_tmap(
        input_file: Path,
        output_dir: Path,
        evaluation_file: Path = None,
) -> None:
    if evaluation_file:
        counts_dict = read_counts(evaluation_file)
        counts_list = []
    lines = read_txt(input_file)[1:]
    smiles, labels = [], []
    for line in lines:
        *bbs, s = line.split('\t')[1:-1]
        smiles += [s]
        label = '-'.join(bbs)
        labels += [label]
        if evaluation_file:
            counts_list += [counts_dict[label]]
    os.chdir(output_dir)
    compute_tree(smiles, labels, counts_list)


if __name__ == '__main__':

    input_file = 'data/smiles_nf.txt'
    output_dir = Path('tmap')
    evaluation_file = 'data/counts.txt'

    run_tmap(input_file, output_dir, evaluation_file)

