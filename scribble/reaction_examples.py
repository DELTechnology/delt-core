from rdkit import Chem
from rdkit.Chem import AllChem

# CuAAC reaction SMARTS (alkyne + azide → 1,2,3‑triazole)
cuaac = AllChem.ReactionFromSmarts(
    "[CX2:1]#[CX2:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N:3][N:4]=[N:5]1")

# Suzuki–Miyaura coupling (aryl iodide + aryl boronic acid → biaryl)
suzuki = AllChem.ReactionFromSmarts(
    "[cX3:1][I].[#6:2][BX3]>>[cX3:1][#6:2]")

# starting molecules
smiles_1 = "OC(=O)C1=CC(=CN=C1)C#C"
smiles_2 = "Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1"
smiles_3 = "BrC1=NC=C(OCC#C)C=C1"
smiles_4 = "OB(O)c1cc(ccc1Cl)C#N"

mol1 = Chem.MolFromSmiles(smiles_1)
mol2 = Chem.MolFromSmiles(smiles_2)
mol3 = Chem.MolFromSmiles(smiles_3)
mol4 = Chem.MolFromSmiles(smiles_4)

# apply CuAAC to two different alkyne partners
product_1 = cuaac.RunReactants((mol1, mol2))[0][0]
product_2 = cuaac.RunReactants((mol3, mol2))[0][0]

# apply Suzuki–Miyaura coupling between the aryl iodide in each CuAAC product and the boronic acid
product_3_from_1 = suzuki.RunReactants((product_1, mol4))[0][0]
product_3_from_2 = suzuki.RunReactants((product_2, mol4))[0][0]

# convert products to SMILES
print("product_1:", Chem.MolToSmiles(product_1))
print("product_2:", Chem.MolToSmiles(product_2))
print("product_3_from_product_1:", Chem.MolToSmiles(product_3_from_1))
print("product_3_from_product_2:", Chem.MolToSmiles(product_3_from_2))
