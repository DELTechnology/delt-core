from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.AllChem import GetRDKitFPGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from sklearn.manifold import TSNE

from utils import read_txt


def read_counts(
        file: Path,
) -> dict:
    counts = defaultdict(int)
    lines = read_txt(file)[1:]
    for line in lines:
        count, bb1, bb2 = line.strip().split('\t')
        counts[f'{bb1}-{bb2}'] = int(count)
    return counts


def compute_fingerprint(
        mol: Chem.rdchem.Mol,
        num_bits: int = 2048,
) -> ExplicitBitVect:
    fpgen = GetRDKitFPGenerator(fpSize=num_bits)
    return fpgen.GetFingerprint(mol)


def compute_similarity(
        bv1: ExplicitBitVect,
        bv2: ExplicitBitVect,
) -> float:
    return DataStructs.TanimotoSimilarity(bv1, bv2)


def plot_fingerprints(
        fingerprints: list,
        counts: list,
        output_file: str = 'fingerprints.png',
) -> None:
    cmap = plt.cm.YlGnBu
    plt.figure(figsize=(10, 5))
    x, y = zip(*fingerprints)
    sc = plt.scatter(x, y, c=counts, edgecolor='k', s=20.0, alpha=1.0, cmap=cmap)
    plt.colorbar(sc, label='Count')
    plt.savefig(output_file, dpi=300)


def run_tsne(
        data: list,
        counts: list,
) -> None:
    smiles = [s.split('\t')[-2] for s in data]
    mols = [Chem.MolFromSmiles(m) for m in smiles]
    fps = [compute_fingerprint(mol) for mol in mols]
    tsne = TSNE(n_components=2, random_state=42)
    fps_2d = tsne.fit_transform(np.array(fps))
    plot_fingerprints(fps_2d, counts)



if __name__ == '__main__':

    input_file = 'smiles_nf.txt'
    counts_file = 'counts.tsv'
    
    num_bb1, num_bb2 = (3, 548)
    data = read_txt(input_file)[1:]
    data = data[:num_bb1 * num_bb2]
    counts = read_counts(counts_file)
    labels = [f'{bb1}-{bb2}' for bb1 in range(1, num_bb1 + 1) for bb2 in range(1, num_bb2 + 1)]
    counts_list = [counts[label] for label in labels]

    run_tsne(data, counts_list)

