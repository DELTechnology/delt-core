from rdkit import Chem
from rdkit.Chem import rdChemReactions

SMIRKS = {
    "CuAAC": "[CX2:1]#[CX2;H1:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N:3][N:4]=[N:5]1",
    "CuAACv2": "[CX2:1]#[CX2:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N:3][N:4]=[N:5]1",
    "SR":    "[#6:1][$([NX2-][NX2+]#[NX1]),$([NX2]=[NX2+]=[NX1-])]>>[#6:1][N;H2]",
    "SRv2":    "[#6:1][N:2]=[N+:3]=[N-:4]>>[#6:1][N;H2]",
    "ABF":   "[CX3:1](=[O:2])[OX2;H1].[N;H2:3]>>[CX3:1](=[O:2])[N:3]",
}

smi_a = "N"
# smi_b = 'Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1'
# smi_b = '[N-]=[N+]=NC(C(O)=O)Cc1cc(I)ccc1'
smi_b = '[N-]=[N+]=NC(C(O)=O)Cc1c(I)cccc1'
rxn_name = 'SR'
rxn = rdChemReactions.ReactionFromSmarts(SMIRKS[rxn_name])
mols = tuple(m for m in (Chem.MolFromSmiles(smi_a), Chem.MolFromSmiles(smi_b)) if m)
products = rxn.RunReactants([mols[1]])
pmol = products[0][0]
Chem.MolToSmiles(pmol, canonical=True)

def run_reaction(rxn_name, smi_a, smi_b=None):
    rxn = rdChemReactions.ReactionFromSmarts(SMIRKS[rxn_name])
    mols = tuple(m for m in (Chem.MolFromSmiles(smi_a), Chem.MolFromSmiles(smi_b)) if m)
    prods = set()
    for tup in rxn.RunReactants(mols):
        for pmol in tup:
            try: Chem.SanitizeMol(pmol, catchErrors=True)
            except: pass
            prods.add(Chem.MolToSmiles(pmol, canonical=True))
    return sorted(prods)