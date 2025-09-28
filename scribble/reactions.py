import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from delt_core.cli.library.api import get_reaction_graph, visualize_reaction_graph, find_next_reaction, \
    perform_reaction, complete_reaction_graph

compounds = {
    's1_1': dict(smiles='OC(=O)C1=CC(=CN=C1)C#C'),
    's2_1': dict(smiles='OB(O)c1cc(ccc1Cl)C#N'),
    '1': dict(smiles='Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1'),
    '2': dict(smiles='[N-]=[N+]=NC(C(O)=O)Cc1cc(I)ccc1'),
    '3': dict(smiles='[N-]=[N+]=NC(C(O)=O)Cc1c(I)cccc1'),
    '4': dict(smiles='NC(Cc1ccc(I)cc1)C(=O)O'),
    '5': dict(smiles='NC(Cc1cccc(I)c1)C(=O)O'),
    '6': dict(smiles='NC(Cc1ccccc1I)C(=O)O'),
}

reactions = {
    'CuAAC': dict(smirks='[CX2:1]#[CX2;H1:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N-0:3][N-0:4]=[N-0:5]1'),
    'SR': dict(smirks='[#6:1][$([NX2-][NX2+]#[NX1]),$([NX2]=[NX2+]=[NX1-])]>>[#6:1][N;H2]'),
    'ABF': dict(smirks='[CX3:1](=[O:2])[OX2;H1].[N;H2:4]>>[CX3:1](=[O:2])[N;H:4]'),
    'Suz': dict(smirks='[cX3:1][I].[#6:2][BX3]>>[cX3:1][#6:2]'),
}

steps = [
    ('s1_1', 'CuAAC'),
    ('1', 'CuAAC'),
    ('CuAAC', 'product_1'),
    ('product_1', 'Suz'),
    ('s2_1', 'Suz'),
    ('Suz', 'product_2'),
]

products = {
    'product_1': dict(smiles=None),
    'product_2': dict(smiles=None),
}

G = get_reaction_graph(steps=steps, reactions=reactions, compounds=compounds, products=products)
ax = visualize_reaction_graph(G=G)
ax.figure.show()

G = complete_reaction_graph(G=G)

next_reaction = find_next_reaction(G)
reactants = [compounds.get(r).get('smiles') for r in next_reaction['reactants']]
smirks = reactions.get(next_reaction['reaction']).get('smirks')

mols = [Chem.MolFromSmiles(i) for i in reactants]
rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=True)
rxn = rdChemReactions.ReactionFromSmarts(smirks)

products = perform_reaction(smirks=smirks, reactants=reactants)

nf = pd.read_excel('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/libraries/NF.xlsx', sheet_name='step1')
cuaac = nf[nf.ReactionType == 'CuAAC']
list_of_reactions = cuaac[['SMILES', 'ScaffoldID', 'ReactionType']].to_dict('records')
for i, reaction in enumerate(list_of_reactions):
    rt = reaction['ReactionType']
    smirks = reactions[rt]['smirks']
    scaffold_id = reaction['ScaffoldID']

    r1 = compounds[str(scaffold_id)]['smiles']
    r2 = reaction['SMILES']

    reactants = [r2, r1]
    products = perform_reaction(smirks=smirks, reactants=reactants)
    if len(products) == 0:
        reactants = [r1, r2]
        products = perform_reaction(smirks=smirks, reactants=reactants)

    if len(products) == 0:
        print(f'No product for reaction: {i}')


nf = pd.read_excel('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/libraries/NF.xlsx', sheet_name='step1')
sr = nf[(nf.ReactionType.str.startswith('SR')) & (nf.ScaffoldID == 1)]
sr = nf[(nf.ReactionType.str.startswith('SR')) & (nf.ScaffoldID == 2)]
sr = nf[(nf.ReactionType.str.startswith('SR')) & (nf.ScaffoldID == 3)]
list_of_reactions = sr[['SMILES', 'ScaffoldID', 'ReactionType']].to_dict('records')
ps = set()
for i, reaction in enumerate(list_of_reactions):
    rt = 'SR'
    smirks = reactions[rt]['smirks']
    scaffold_id = reaction['ScaffoldID']

    r1 = compounds[str(scaffold_id)]['smiles']
    r2 = reaction['SMILES']

    reactants = [r1]
    products = perform_reaction(smirks=smirks, reactants=reactants)
    if len(products) == 0:
        print(f'No product for reaction: {i}')
    ps.update(products)

nf = pd.read_excel('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/libraries/NF-v2.xlsx', sheet_name='step1')
sf4 = nf[nf.ScaffoldID == 4]
list_of_reactions = sf4[['SMILES', 'ScaffoldID', 'ReactionType']].to_dict('records')
for i, reaction in enumerate(list_of_reactions):
    rt = reaction['ReactionType']
    smirks = reactions[rt]['smirks']
    scaffold_id = reaction['ScaffoldID']

    r1 = compounds[str(scaffold_id)]['smiles']
    r2 = reaction['SMILES']

    reactants = [r2, r1]
    products = perform_reaction(smirks=smirks, reactants=reactants)
    if len(products) == 0:
        reactants = [r1, r2]
        products = perform_reaction(smirks=smirks, reactants=reactants)

    if len(products) == 0:
        print(f'No product for reaction: {i}')

nf = pd.read_excel('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/libraries/NF-v2.xlsx', sheet_name='step1')
sf5 = nf[nf.ScaffoldID == 5]
list_of_reactions = sf5[['SMILES', 'ScaffoldID', 'ReactionType']].to_dict('records')
for i, reaction in enumerate(list_of_reactions):
    rt = reaction['ReactionType']
    smirks = reactions[rt]['smirks']
    scaffold_id = reaction['ScaffoldID']

    r1 = compounds[str(scaffold_id)]['smiles']
    r2 = reaction['SMILES']

    reactants = [r2, r1]
    products = perform_reaction(smirks=smirks, reactants=reactants)
    if len(products) == 0:
        reactants = [r1, r2]
        products = perform_reaction(smirks=smirks, reactants=reactants)

    if len(products) == 0:
        print(f'No product for reaction: {i}')

nf = pd.read_excel('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/DECLT-DB/libraries/NF-v2.xlsx', sheet_name='step1')
sf6 = nf[nf.ScaffoldID == 6]
list_of_reactions = sf6[['SMILES', 'ScaffoldID', 'ReactionType']].to_dict('records')
for i, reaction in enumerate(list_of_reactions):
    rt = reaction['ReactionType']
    smirks = reactions[rt]['smirks']
    scaffold_id = reaction['ScaffoldID']

    r1 = compounds[str(scaffold_id)]['smiles']
    r2 = reaction['SMILES']

    reactants = [r2, r1]
    products = perform_reaction(smirks=smirks, reactants=reactants)
    if len(products) == 0:
        reactants = [r1, r2]
        products = perform_reaction(smirks=smirks, reactants=reactants)

    if len(products) == 0:
        print(f'No product for reaction: {i}')
