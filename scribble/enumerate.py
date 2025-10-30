import networkx as nx
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions

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


def visualize_reaction_graph(G: nx.DiGraph) -> plt.Axes:
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    compounds = [n for n, d in G.nodes(data=True) if d.get("type") == "compound"]
    reactions = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    products = [n for n, d in G.nodes(data=True) if d.get("type") == "product"]

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")

    # Draw compounds (blue circles)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=compounds,
        node_color="lightblue",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw compounds (blue circles)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=products,
        node_color="salmon",
        node_shape="o",  # circle
        node_size=500,
        ax=ax
    )

    # Draw reactions (green squares)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=reactions,
        node_color="lightgreen",
        node_shape="s",  # square
        node_size=600,
        ax=ax
    )

    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True)

    fig.show()
    return ax


def find_next_reaction(G: nx.DiGraph):
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == "reaction"]
    for node in reaction_nodes:
        preds = sorted(G.predecessors(node))
        succ, = sorted(G.successors(node))  # note: currently only one product per reaction

        if G.nodes[succ]['smiles'] is None and all([G.nodes[i]['smiles'] is not None for i in preds]):
            return {'reactants': preds, 'reaction': node, 'product': succ}

    return None


def perform_reaction(smirks: str, reactants: list[str]) -> list[str]:
    mols = [Chem.MolFromSmiles(i) for i in reactants]
    rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=True)

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


def complete_reaction_graph(G: nx.DiGraph) -> nx.DiGraph:
    while True:

        try:
            next_reaction = find_next_reaction(G)
            if next_reaction is None:
                break

            smirks = G.nodes[next_reaction['reaction']]['smirks']
            reactants = [G.nodes[i]['smiles'] for i in next_reaction['reactants']]
            products = perform_reaction(smirks, reactants)

            # assert len(products) == 1, f"Expected one product, got {products}"
            product = {next_reaction['product']: dict(smiles=products[0])}
            nx.set_node_attributes(G, product)
        except Exception as e:
            print(f"Error processing reaction {next_reaction}: {e}")
            break

    return G


compounds = {
    'alkyne': dict(smiles="CC#C"),
    'azide': dict(smiles="N=N#Nc1ccccc1"),
}

# DOCS:
# - https://rdkit.blogspot.com/2015/01/chemical-reaction-notes-i.html
reactions = {
    'CuAAC': dict(smirks="[*:7][C:1]#[C:2][*:8].[N:3]=[N+:4]=[N-:5]>>[*:7][c:1]1[n:3][n:4][n:5][c:2]1[*:8]"),
    'CuAAC-terminal': dict(smirks="[*:7][C:1]#[C:2][H].[N:3]=[N+:4]=[N-:5]>>[*:7][c:1]1[n:3][n:4][n:5][c:2]1"),
}

products = {
    'product_1': dict(smiles=None),
}

steps = [
    ('alkyne', 'CuAAC'),
    ('azide', 'CuAAC'),
    ('CuAAC', 'product_1'),
]

G = get_reaction_graph(steps=steps, reactions=reactions, compounds=compounds, products=products)
ax = visualize_reaction_graph(G)

G = complete_reaction_graph(G)
smiles = G.nodes['product_1']['smiles']

mols = [Chem.MolFromSmiles(smiles)]
img = Draw.MolsToGridImage(mols, molsPerRow=1)
plt.imshow(img).figure.show()