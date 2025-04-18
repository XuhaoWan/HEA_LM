#!/usr/bin/env python
# -*- coding:utf-8 -*-
#author: xhwan

import math

radius_dict = {
    'Co': 126,
    'Cr': 127,
    'Fe': 125,
    'Ni': 121,
    'Cu': 138,
    'Al': 118,
    'Ti': 136,
    'Ga': 126,
    'Mn': 139,
    'Mo': 145
}

metal_point_dict = {
    'Co': 1768,
    'Cr': 2180,
    'Fe': 1811,
    'Ni': 1728,
    'Cu': 1357.77,
    'Al': 933.47,
    'Ti': 1941,
    'Ga': 302.91,
    'Mn': 1519,
    'Mo': 2896
}

energy_fcc_dict = {
    'Co': -6.278748125,
    'Cr': -8.933871679,
    'Fe': -7.330305539,
    'Ni': -5.234804971,
    'Cu': -3.887359329,
    'Al': -3.881197387,
    'Ti': -8.12662124,
    'Ga': -2.92122599,
    'Mn': -8.274325915,
    'Mo': -11.02389704
}

energy_bcc_dict = {
    'Co': -6.744816978,
    'Cr': -9.845993166,
    'Fe': -7.878695222,
    'Ni': -5.64281086,
    'Cu': -4.167520352,
    'Al': -3.340218292,
    'Ti': -7.580988719,
    'Ga': -2.155807907,
    'Mn': -8.964668019,
    'Mo': -10.60220741
}


def parse_alloy_composition(composition_str):

    elements = composition_str.replace(',', ' ').split()
    

    composition_dict = {}
    
    for element in elements:

        for i in range(len(element)):
            if element[i].isdigit():
                element_symbol = element[:i]
                e_N = element[i:]
                break
        
        e_compos = int(e_N)
        
        composition_dict[element_symbol] = e_compos
    
    return composition_dict


def calculate_average_radius(composition_str, radius_dict):
    composition = parse_alloy_composition(composition_str)
    total_radius_contribution = 0.0
    total_percentage = 0.0
    
    for element, percentage in composition.items():
        if element in radius_dict:
            total_radius_contribution += percentage * radius_dict[element]
            total_percentage += percentage
        else:
            print(f"element {element} not found.")
    
    if total_percentage == 0:
        return None
    
    return total_radius_contribution / total_percentage, total_percentage


def calculate_average_metal_point(composition_str, metal_point_dict):
    composition = parse_alloy_composition(composition_str)
    total_TM_contribution = 0.0
    total_percentage = 0.0
    
    for element, percentage in composition.items():
        if element in metal_point_dict:
            total_TM_contribution += percentage * metal_point_dict[element]
            total_percentage += percentage
        else:
            print(f"element {element} not found.")
    
    if total_percentage == 0:
        return None
    return total_TM_contribution / total_percentage, total_percentage


def calculate_radius_delta(composition_str, radius_dict):
    composition = parse_alloy_composition(composition_str)
    Rave, e_N = calculate_average_radius(composition_str, radius_dict)
    total_element_contribution = 0.0

    
    for element, percentage in composition.items():
        if element in radius_dict:
            total_element_contribution += percentage * math.pow((1-radius_dict[element]/Rave), 2)
            #print(element, math.pow((1-radius_dict[element]/Rave), 2))
        else:
            print(f"element {element} not found.")

    delta = math.sqrt(total_element_contribution/e_N)

    return delta


def calculate_sigma(lat_type, energy, composition_str, energy_fcc_dict, energy_bcc_dict):
    lat = lat_type
    composition = parse_alloy_composition(composition_str)
    TM, e_N = calculate_average_metal_point(composition_str, metal_point_dict)
    Epure_all = read(raw.csv)

    Epure = Epure_all/e_N
    dH = energy/e_N - Epure
    dS = -ra * (8.314*6.24151/602214)
    sigma = TM*dS/abs(dH)

    return sigma, dS, dH


