#!/usr/bin/env python
# -*- coding:utf-8 -*-
#author: xhwan

import numpy as np
import pandas as pd
import networkx as nx
from ase.io import read
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt


fcc_lattice_constants = {
    "Co": 3.544, "Ga": "u", "Ni": 3.524, "Cu": 4.05, "Ti": "u",
    "Cr": "u", "Fe": 3.65, "Al": 4.05, "Mn": 3.645, "Mo": "u"
}

bcc_lattice_constants = {
    "Co": 2.821, "Ga": "u", "Ni": "u", "Cu": "u", "Ti": 3.31,
    "Cr": 2.88, "Fe": 2.87, "Al": "u", "Mn": 2.89, "Mo": 3.15
}

ELEMENT_DATA = {
    'Element': ['Co', 'Cr', 'Fe', 'Ni', 'Cu', 'Al', 'Ti', 'Ga', 'Mn', 'Mo'],
    'AtomicNumber': [27, 24, 26, 28, 29, 13, 22, 31, 25, 42],
    'AtomicRadius': [126, 127, 125, 121, 138, 118, 136, 126, 139, 145],
    'AtomicMass': [58.93, 51.99, 55.84, 58.69, 63.55, 26.98, 47.87, 69.72, 54.94, 95.96],
    'Electronegativity': [1.88, 1.66, 1.83, 1.91, 1.90, 1.61, 1.54, 1.81, 1.55, 2.16],
    'MeltingPoint': [1768, 2180, 1811, 1728, 1357, 933, 1941, 303, 1519, 2896],
    'BoilingPoint': [3200, 2944, 3134, 3186, 2835, 2470, 3560, 2477, 2334, 4912],
    'Period': [4, 4, 4, 4, 4, 3, 4, 4, 4, 5],
    'Group': [9, 6, 8, 10, 11, 13, 4, 13, 7, 6],
    'OxidationMax': [5, 6, 7, 4, 4, 3, 4, 3, 7, 6],
    'OxidationMin': [-3, -4, -4, -2, -2, -2, -2, -5, -3, -4],
    'IonizationEnergy': [7.881, 6.767, 7.903, 7.640, 7.727, 5.985, 6.828, 5.999, 7.434, 7.092],
    'Electronaffinity': [0.660, 0.666, 0.163, 1.160, 1.227, 0.440, 0.079, 0.300, 0, 0.745]
}

df_elements = pd.DataFrame(ELEMENT_DATA)
df_elements.set_index('Element', inplace=True)


def encode_onehot(value, categories):
    vec = np.zeros(len(categories))
    vec[categories.index(value)] = 1
    return vec


def calculate_lattice_mismatch(element1, element2, total_atoms):
    
    if total_atoms == 128:
        lattice_constants = bcc_lattice_constants
    elif total_atoms == 108:
        lattice_constants = fcc_lattice_constants
    else:
        raise ValueError("total atom number must be 128 (BCC) or 108 (FCC)")

    a1 = lattice_constants.get(element1)
    a2 = lattice_constants.get(element2)

    if a1 == "u" or a2 == "u":
        return 0.0

    mismatch = abs(a1 - a2) / ((a1 + a2) / 2)
    return mismatch


def get_atomic_features(element):
    data = df_elements.loc[element]

    onehot_atomic_number = encode_onehot(data['AtomicNumber'], list(range(1, 43)))
    onehot_period = encode_onehot(data['Period'], list(range(1, 8)))
    onehot_group = encode_onehot(data['Group'], list(range(1, 19)))

    continuous_features = data[['AtomicRadius', 'AtomicMass', 'Electronegativity',
                                'MeltingPoint', 'BoilingPoint', 'IonizationEnergy', 
                                'Electronaffinity']].values
    oxidation_states = [data['OxidationMax'], data['OxidationMin']]

    return np.concatenate([onehot_atomic_number, onehot_period, onehot_group,
                           continuous_features, oxidation_states])


def build_graph_from_cif(cif_file, cutoff):
    atoms = read(cif_file)
    positions = atoms.get_positions()
    atom_types = atoms.get_chemical_symbols()
    total_atoms = len(atoms)
    cell_volume = atoms.get_volume()

    total_mass = np.sum(atoms.get_masses())    
    density = total_mass / cell_volume

    G = nx.Graph(volume=cell_volume, density=density)

    for i, (atom, pos) in enumerate(zip(atom_types, positions)):
        features = get_atomic_features(atom)
        G.add_node(i, features=features)

    distances = atoms.get_all_distances(mic=True)

    for i in range(len(atom_types)):
        for j in range(i + 1, len(atom_types)):
            if distances[i, j] < cutoff:
                electronegativity_diff = abs(
                    df_elements.loc[atom_types[i], 'Electronegativity'] -
                    df_elements.loc[atom_types[j], 'Electronegativity']
                )
                mismatch = calculate_lattice_mismatch(atom_types[i], atom_types[j], total_atoms)
                G.add_edge(i, j, length=distances[i, j], 
                           electronegativity_diff=electronegativity_diff, mismatch=mismatch)

    for node in G.nodes:
        G.nodes[node]['features'] = np.append(G.nodes[node]['features'], G.degree[node])


    return G


def graph_to_vector(G, max_nodes=128, max_edges=896):
    vector = []

    vector.append(G.graph.get('volume', 0))
    vector.append(G.graph.get('density', 0))

    node_features = []
    for i in range(max_nodes):
        if i in G:
            node_features.extend(G.nodes[i].get('features', np.zeros(77)))  
        else:
            node_features.extend(np.zeros(77))

    vector.extend(node_features)

    edge_features = []
    edges = list(G.edges(data=True))
    for i in range(max_edges):
        if i < len(edges):
            edge_data = edges[i][2]
            edge_features.extend([
                edge_data.get('length', 0),
                edge_data.get('electronegativity_diff', 0),
                edge_data.get('mismatch', 0)
            ])
        else:
            edge_features.extend([0, 0, 0])

    vector.extend(edge_features)
    return np.array(vector)


if __name__ == "__main__":
    cif_file = 'hea.cif'

    # 构建图
    graph  = build_graph_from_cif(cif_file, cutoff=3.5)

#    connected_edges = graph.edges(0)
#    for edge in connected_edges:
#        edge_data = graph[edge[0]][edge[1]] 
#        print(f"Edge {edge}: Data = {edge_data}")
#    print(graph.nodes[0].get('features'))
#    print(graph.graph)
#    edge_features = []
#    edges = list(graph.edges(data=True))
#    print(edges)
#    nx.draw(graph)
#    plt.show()
    vec = graph_to_vector(graph)
    df = pd.DataFrame([vec])
