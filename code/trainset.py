#!/usr/bin/env python
# -*- coding:utf-8 -*-
#author: xhwan

from sigma import parse_hea_phase
from mkhea import generate_random_structure, save_cif, parse_input, generate_alloy
import pandas as pd
import os



def generate_all_heas(num_fcc=600, num_bcc=600):
    os.makedirs("trainheas", exist_ok=True)

    records = []

    for i in range(1, num_fcc + num_bcc + 1):
        lattice_type = 'fcc' if i <= num_fcc else 'bcc'
        
        elements = generate_alloy(lattice_type, 4, 6)
        print(i, elements)
        
        structure = generate_random_structure(parse_input(elements))
        filename = os.path.join("trainheas", f"hea{i}.cif")
        save_cif(structure, filename=filename)

        
        includes_ga = 1 if 'Ga' in elements else 0

        record = {
            "HEA_ID": i,
            "Elements": elements,
            "Lattice_Type": lattice_type,
            "Phase": phase,
            "Energy (E)": E,
            "Includes_Ga": includes_ga,
            "Delta": delta,
            "Sigma": sigma
        }
        records.append(record)

    df = pd.DataFrame(records)
    df.to_csv("hea_summary.csv", index=False)

    print("Generation complete. Saved to 'hea_summary.csv'.")

generate_all_heas()
