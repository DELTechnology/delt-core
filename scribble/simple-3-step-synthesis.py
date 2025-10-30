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

    # Draw products (salmon circles)
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

            # Take the first product if multiple are formed
            product = {next_reaction['product']: dict(smiles=products[0])}
            nx.set_node_attributes(G, product)
            print(f"Reaction {next_reaction['reaction']}: {reactants} -> {products[0]}")
        except Exception as e:
            print(f"Error processing reaction {next_reaction}: {e}")
            break

    return G


# Simple three-step reaction: alkene -> alcohol -> ketone -> alcohol
compounds = {
    'ethene': dict(smiles="C=C"),  # Starting alkene
    'water': dict(smiles="O"),  # Water for hydration
    'oxidant': dict(smiles="[O]"),  # Atomic oxygen for oxidation
    'reductant': dict(smiles="[H][H]"),  # Hydrogen gas for reduction
}

reactions = {
    # Step 1: Alkene hydration - ethene + water -> ethanol
    'hydration': dict(smirks="[C:1]=[C:2].[O:3]>>[C:1][C:2][O:3]"),

    # Step 2: Alcohol oxidation - ethanol -> acetaldehyde
    'oxidation': dict(smirks="[C:1][C:2][O:3].[O:4]>>[C:1][C:2]=[O:4].[O:3]"),

    # Step 3: Aldehyde reduction - acetaldehyde + H2 -> ethanol
    # One H goes to carbon, one H goes to oxygen to form C-OH
    'reduction': dict(smirks="[C:1][C:2]=[O:3].[H:4][H:5]>>[C:1][C:2]([H:4])[O:3][H:5]"),
}

products = {
    'intermediate_1': dict(smiles=None),  # Ethanol after hydration
    'intermediate_2': dict(smiles=None),  # Acetaldehyde after oxidation
    'final_product': dict(smiles=None),  # Ethanol after reduction
}

steps = [
    # Step 1: Hydration
    ('ethene', 'hydration'),
    ('water', 'hydration'),
    ('hydration', 'intermediate_1'),

    # Step 2: Oxidation
    ('intermediate_1', 'oxidation'),
    ('oxidant', 'oxidation'),
    ('oxidation', 'intermediate_2'),

    # Step 3: Reduction
    ('intermediate_2', 'reduction'),
    ('reductant', 'reduction'),
    ('reduction', 'final_product'),
]

# Create and visualize the reaction graph
G = get_reaction_graph(steps=steps, reactions=reactions, compounds=compounds, products=products)
ax = visualize_reaction_graph(G)

# Complete the reactions
G = complete_reaction_graph(G)

# Display final products
print("\nFinal products:")
for node, data in G.nodes(data=True):
    if data.get('type') in ['product'] and data.get('smiles'):
        print(f"{node}: {data['smiles']}")

# Visualize the molecular structures of all products
product_smiles = []
product_names = []
for node, data in G.nodes(data=True):
    if data.get('type') == 'product' and data.get('smiles'):
        product_smiles.append(data['smiles'])
        product_names.append(node)

if product_smiles:
    mols = [Chem.MolFromSmiles(smiles) for smiles in product_smiles]
    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=product_names)
    plt.figure(figsize=(10, 6))
    plt.imshow(img)
    plt.axis('off')
    plt.title('Product Structures')
    plt.show()