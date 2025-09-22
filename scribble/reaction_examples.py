from rdkit import Chem
from rdkit.Chem import AllChem
from delt_core.assembly.reaction_graph import visualize_reaction_graph
import networkx as nx

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
}

reactions = {
    'CuAAC':  dict(smirks="[CX2:1]#[CX2:2].[N:3]=[N+:4]=[N-:5]>>[C:1]1=[C:2][N:3][N:4]=[N:5]1"),
    'Suzuki': dict(smirks="[cX3:1][I].[#6:2][BX3]>>[cX3:1][#6:2]")
}

products = {'product_1': dict(smiles=None),
            'product_2': dict(smiles=None)}

steps = [
    ('smiles_1', 'CuAAC'),
    ('smiles_2', 'CuAAC'),
    ('CuAAC', 'product_1'),
    ('product_1', 'Suzuki'),
    ('smiles_4', 'Suzuki'),
    ('Suzuki', 'product_2'),
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
            # sanitize and canonicalize
            Chem.SanitizeMol(pmol, catchErrors=True)
            products.add(Chem.MolToSmiles(pmol, canonical=True))

    return sorted(products)

products = perform_reaction(smirks, reactants)



# mol1 = Chem.MolFromSmiles(compounds['smiles_1'])
# mol2 = Chem.MolFromSmiles(compounds['smiles_2'])
# mol3 = Chem.MolFromSmiles(compounds['smiles_3'])
# mol4 = Chem.MolFromSmiles(compounds['smiles_4'])
#
# # apply CuAAC to two different alkyne partners
# product_1 = cuaac.RunReactants((mol1, mol2))[0][0]
# product_2 = cuaac.RunReactants((mol3, mol2))[0][0]
#
# # apply Suzuki–Miyaura coupling between the aryl iodide in each CuAAC product and the boronic acid
# product_3_from_1 = suzuki.RunReactants((product_1, mol4))[0][0]
# product_3_from_2 = suzuki.RunReactants((product_2, mol4))[0][0]
#
# # convert products to SMILES
# print("product_1:", Chem.MolToSmiles(product_1))
# print("product_2:", Chem.MolToSmiles(product_2))
# print("product_3_from_product_1:", Chem.MolToSmiles(product_3_from_1))
# print("product_3_from_product_2:", Chem.MolToSmiles(product_3_from_2))
