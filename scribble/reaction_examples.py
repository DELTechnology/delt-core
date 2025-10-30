import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem

from delt_core.assembly.reaction_graph import visualize_reaction_graph


def get_reaction_graph(steps: list, reactions: dict, compounds: dict, products: dict) -> nx.DiGraph:
    G = nx.DiGraph()

    G.add_edges_from(steps)

    attrs = {k: {**v, 'type': 'reaction'} for k, v in reactions.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'compound'} for k, v in compounds.items()}
    nx.set_node_attributes(G, attrs)

    attrs = {k: {**v, 'type': 'product'} for k, v in products.items()}
    nx.set_node_attributes(G, attrs)

    return G


compounds = {
    'smiles_1': dict(smiles="OC(=O)C1=CC(=CN=C1)C#C"),
    'smiles_2': dict(smiles="Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1"),
    'smiles_3': dict(smiles="BrC1=NC=C(OCC#C)C=C1"),
    'smiles_4': dict(smiles="OB(O)c1cc(ccc1Cl)C#N"),
    'smiles_5': dict(smiles="c1ccccc1"),
    'alkyne': dict(smiles="CC#C"),  # propyne (R=CH3, terminal H present)
    'azide': dict(smiles="c1ccccc1N=[N+]=[N-]"), # phenyl azide
    'smiles_a': dict(smiles='C(=O)O'),
    'smiles_b': dict(smiles="CNC"),
}

# DOCS:
# - https://rdkit.blogspot.com/2015/01/chemical-reaction-notes-i.html
reactions = {
    'CuAAC': dict(smirks="[CX2:1]#[CX2:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N:3][N:4]=[N:5]1"),
    'Suzuki': dict(smirks="[cX3:1][I].[#6:2][BX3]>>[cX3:1][#6:2]"),
    'Suzuki-v2': dict(smirks="(O)(O)([#6:1]).[Cl,Br,I,$(OS(=O)(=O)C(F)(F)F)][#6:2]>>[#6:1]-[#6:2]"),
    'CuAAC-v2': dict(
        smirks="[#6,#7,#8,#16,*:1]-[C:2]#[C:3]-[H:4].[N:5]=[N:6]=[N:7]-[*:8]>>[n:5]1[c:3][n:6][n:7][c:2]1-[*:8].[#6,#7,#8,#16,*:1]>>[n:5]1[c:3][n:6][n:7][c:2]1-[*:8]-[*:1]"),
    'CuAAC-v3': dict(smirks="[*:1]-[C:2]#[C:3]-[H:4].[N:5]=[N:6]=[N:7]-[*:8]>>[n:5]1[c:3][n:6][n:7][c:2]1-[*:8]-[*:1]"),
    'CuAAC-v4': dict(smirks="[N-:1]=[N+:2]=[N:3]-[*:4]>>[*:4]-[N:3]=[N+:2]=[N-:1]"),
    'CuAAC-latest': dict(smirks="[C:1]#[C:2].[N:3]=[N+:4]=[N-:5]>>[c:1]1[n:3][n:4][n:5][c:2]1"),
    'tutorial-rxn': dict(smirks='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
}

products = {
    'product_1': dict(smiles=None),
    'product_2': dict(smiles=None),
    'product_v2': dict(smiles=None),
    'product_v3': dict(smiles=None),
    'product_v4': dict(smiles=None),
    'product_c': dict(smiles=None),
}

steps = [
    ('smiles_1', 'CuAAC'),
    ('smiles_2', 'CuAAC'),
    ('CuAAC', 'product_1'),
    ('product_1', 'Suzuki'),
    ('smiles_4', 'Suzuki'),
    ('Suzuki', 'product_2'),
]

steps = [
    ('alkyne', 'CuAAC-v2'),
    ('azide', 'CuAAC-v2'),
    ('CuAAC-v2', 'product_v2'),
]

steps = [
    ('alkyne', 'CuAAC-v3'),
    ('azide', 'CuAAC-v3'),
    ('CuAAC-v3', 'product_v3'),
]

steps = [
    ('alkyne', 'CuAAC-latest'),
    ('azide', 'CuAAC-latest'),
    ('CuAAC-latest', 'product_v4'),
]

steps = [
    ('smiles_a', 'tutorial-rxn'),
    ('smiles_b', 'tutorial-rxn'),
    ('tutorial-rxn', 'product_c'),
]

G = get_reaction_graph(steps=steps, reactions=reactions, compounds=compounds, products=products)
ax = visualize_reaction_graph(G)


def find_next_reaction(G: nx.DiGraph):
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    for node in reaction_nodes:
        preds = sorted(G.predecessors(node))
        succs, = sorted(G.successors(node))  # note: currently only one product per reaction

        if G.nodes[succs]['smiles'] is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
            return {'reactants': preds, 'reaction': node, 'product': succs}

    return None


next_reaction = find_next_reaction(G)
smirks = reactions[next_reaction['reaction']]['smirks']
reactants = [compounds[i]['smiles'] for i in next_reaction['reactants']]

def perform_reaction(smirks: str, reactants: list[str]) -> list[str]:
    mols = [Chem.MolFromSmiles(i) for i in reactants]
    rxn = AllChem.ReactionFromSmarts(smirks)

    # note: order matters
    product_sets = rxn.RunReactants([*mols])

    products = set()
    for tup in product_sets:
        for pmol in tup:
            # break
            # sanitize and canonicalize
            # flags = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            # Chem.SanitizeMol(m, sanitizeOps=flags)
            # Chem.SetAromaticity(m)

            Chem.SanitizeMol(pmol)
            products.add(Chem.MolToSmiles(pmol, canonical=True, kekuleSmiles=False, isomericSmiles=False))

    return sorted(products)

prods = perform_reaction(smirks, reactants)
product = {next_reaction['product']: dict(smiles=prods[0])}
nx.set_node_attributes(G, product)

def reset_and_aromatize(m):
    m = Chem.Mol(m)
    for a in m.GetAtoms():
        a.SetIsAromatic(False)
    for b in m.GetBonds():
        b.SetIsAromatic(False)
    Chem.SetAromaticity(m)
    return m

def complete_reaction_graph(G: nx.DiGraph) -> nx.DiGraph:
    while True:

        try:
            next_reaction = find_next_reaction(G)
            if next_reaction is None:
                break

            smirks = reactions[next_reaction['reaction']]['smirks']
            reactants = [G.nodes[i]['smiles'] for i in next_reaction['reactants']]
            products = perform_reaction(smirks, reactants)

            assert len(products) == 1, f"Expected one product, got {products}"
            product = {next_reaction['product']: dict(smiles=products[0])}
            nx.set_node_attributes(G, product)
        except Exception as e:
            print(f"Error processing reaction {next_reaction}: {e}")
            break

    return G

m = reset_and_aromatize(pmol)

G = complete_reaction_graph(G)

# from rdkit.Chem import Draw
# from matplotlib import pyplot as plt
# mols = [Chem.MolFromSmiles(i) for i in reactants]
# img = Draw.MolsToGridImage( mols, molsPerRow=6 )
# plt.imshow(img).figure.show()
