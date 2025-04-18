#!/usr/bin/env python
# -*- coding:utf-8 -*-
#author: xhwan


import random
from ase import Atoms
from ase.io import write
from ase.build import bulk, make_supercell



elements = ['Co', 'Cr', 'Fe', 'Ni', 'Cu', 'Al', 'Ti', 'Ga', 'Mn', 'Mo']

fcc_lattice_constants = {
    "Co": 3.544, "Ga": "u", "Ni": 3.524, "Cu": 4.05, "Ti": "u",
    "Cr": "u", "Fe": 3.65, "Al": 4.05, "Mn": 3.645, "Mo": "u"
}

bcc_lattice_constants = {
    "Co": 2.821, "Ga": "u", "Ni": "u", "Cu": "u", "Ti": 3.31,
    "Cr": 2.88, "Fe": 2.87, "Al": "u", "Mn": 2.89, "Mo": 3.15
}



def generate_alloy(lat_type, min_elements=5, max_elements=6):

    if lat_type == 'fcc':
        total_atoms = 108
    elif lat_type == 'bcc':
        total_atoms = 128
    else:
        print('lat type is wrong input')

    num_elements = random.randint(min_elements, max_elements)
    selected_elements = random.sample(elements, num_elements)
    
    min_atoms = int(total_atoms * 0.05)
    max_atoms = int(total_atoms * 0.35)
    
    atom_counts = []
    remaining_atoms = total_atoms
    
    for i in range(len(selected_elements)):
        if i == len(selected_elements) - 1:
            atom_counts.append(remaining_atoms)
        else:
            count = random.randint(min_atoms, min(max_atoms, remaining_atoms - (len(selected_elements) - i - 1) * min_atoms))
            atom_counts.append(count)
            remaining_atoms -= count
    
    alloy_string = ' '.join(f'{element}{count}' for element, count in zip(selected_elements, atom_counts))
    return alloy_string


def parse_input(input_str):
    elements = []
    parts = input_str.split()
    
    for part in parts:
        element = ''.join(filter(str.isalpha, part))
        count = int(''.join(filter(str.isdigit, part)))
        elements.extend([element] * count)
        
    return elements

def get_average_lattice_constant(elements):
    structure_type = None
    if len(elements) == 128:
        structure_type = 'bcc'
    elif len(elements) == 108:
        structure_type = 'fcc'
    else:
        raise ValueError('The input elements number is not satisfied')

    if structure_type == "fcc":
        constants = fcc_lattice_constants
    elif structure_type == "bcc":
        constants = bcc_lattice_constants
    else:
        raise ValueError("Invalid structure type! Choose 'fcc' or 'bcc'.")

    valid_constants = [
        constants[el] for el in elements if constants[el] != "u"
    ]

    if not valid_constants:
        if structure_type == "fcc":
            valid_constants = [3.74383]
        elif structure_type == "bcc":
            valid_constants = [2.98683]
        else:
            raise ValueError("Problem when get lat constant.")
    
    return sum(valid_constants) / len(valid_constants)

def generate_random_structure(elements):
    structure_type = None

    if len(elements) == 128:
        structure_type = 'bcc'
    elif len(elements) == 108:
        structure_type = 'fcc'
    else:
        raise ValueError('The input elements number is not satisfied')

    lattice_constant = get_average_lattice_constant(elements)
    random.shuffle(elements)  

    num_atoms = len(elements)
    
    base_bulk = bulk('X', crystalstructure=structure_type, a=lattice_constant, cubic=True)
    if structure_type == 'bcc':
        hea = make_supercell(base_bulk, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    elif structure_type == 'fcc':
        hea = make_supercell(base_bulk, [[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    else:
        raise ValueError(f"wrong structure_type.")
    hea.set_chemical_symbols(elements[:num_atoms])  

    return hea

def save_cif(structure, filename="hea.cif"):
    write(filename, structure)
    print(f"Structure saved as {filename}")

if __name__ == "__main__":
    user_input = "Ga26 Mn26 Fe26 Ni25 Cu25"  
    elements = parse_input(user_input)

    try:
        alloy_structure = generate_random_structure(elements)
        save_cif(alloy_structure)
    except ValueError as e:
        print(e)
