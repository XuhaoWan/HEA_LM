#!/usr/bin/env python
# -*- coding:utf-8 -*-
#author: xhwan


import os
import pandas as pd
from hea2graph import build_graph_from_cif, graph_to_vector
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


def process_all_heas(cif_folder="trainheas", summary_file="hea_summary.csv", output_file="hea_train_features.csv"):
    hea_summary = pd.read_csv(summary_file)
    
    records = []

    for index, row in hea_summary.iterrows():
        hea_id = row["HEA_ID"]
        phase = row["Phase"]
        lattice_type = row["Lattice_Type"]
        includes_ga = row["Includes_Ga"]
        
        cif_path = os.path.join(cif_folder, f"hea{hea_id}.cif")
        
        if os.path.exists(cif_path):
            graph = build_graph_from_cif(cif_path, cutoff=3.5)
            
            vector = graph_to_vector(graph)
            
            vector = [hea_id, phase, lattice_type, includes_ga] + list(vector)
            
            records.append(vector)
        else:
            print(f"Warning: CIF file {cif_path} not found.")
    
    columns = ["HEA_ID", "Phase", "Lattice_Type", "Includes_Ga"] + [f"Feature_{i+1}" for i in range(len(vector) - 4)]
    df = pd.DataFrame(records, columns=columns)
    df.to_csv(output_file, index=False)

    print(f"Feature extraction complete. Saved to '{output_file}'.")



if __name__ =='__main__':
    data = pd.read_csv('hea_train_features.csv')


    group1 = data[(data['Lattice_Type'] == 'fcc') & (data['Includes_Ga'] == 1)]

    group2 = data[(data['Lattice_Type'] == 'fcc') & (data['Includes_Ga'] == 0)]

    group3 = data[(data['Lattice_Type'] == 'bcc') & (data['Includes_Ga'] == 1)]

    group4 = data[(data['Lattice_Type'] == 'bcc') & (data['Includes_Ga'] == 0)]


    features1 = group1.drop(columns=['HEA_ID', 'Lattice_Type'])
    features2 = group2.drop(columns=['HEA_ID', 'Lattice_Type'])
    features3 = group3.drop(columns=['HEA_ID', 'Lattice_Type'])
    features4 = group4.drop(columns=['HEA_ID', 'Lattice_Type'])


    all_data = pd.concat([features1, features2, features3, features4])

    all_data = all_data.loc[:, (all_data != 0).any(axis=0)]

    scaler = MinMaxScaler()

    all_data = scaler.fit_transform(all_data)

    tsne = TSNE(n_components=2, perplexity = 50, init='pca', learning_rate=50)
    tsne_results = tsne.fit_transform(all_data)

    plt.figure(figsize=(10, 7))
    
    plt.scatter(tsne_results[:len(features1), 0], tsne_results[:len(features1), 1], label='fcc_SS', alpha=1, c='r')
    plt.scatter(tsne_results[len(features1):len(features1)+len(features2), 0], tsne_results[len(features1):len(features1)+len(features2), 1], label='fcc_noSS', c = 'g', alpha=0.5)
    plt.scatter(tsne_results[len(features1)+len(features2):len(features1)+len(features2)+len(features3), 0], tsne_results[len(features1)+len(features2):len(features1)+len(features2)+len(features3), 1], label='bcc_SS', alpha=0.5)
    plt.scatter(tsne_results[len(features1)+len(features2)+len(features3):len(features1)+len(features2)+len(features3)+len(features4), 0], tsne_results[len(features1)+len(features2)+len(features3):len(features1)+len(features2)+len(features3)+len(features4), 1], label='bcc_no_SS', alpha=0.5)
    

            
    plt.legend()
    plt.title("t-SNE of the 4 Groups")
    plt.xlabel("t-SNE Component 1")
    plt.ylabel("t-SNE Component 2")
    plt.savefig('tsne.png')
    plt.show()



    def get_tsne_with_label(group_tsne_results, group_label):
        tsne_df = pd.DataFrame(group_tsne_results, columns=['tSNE_1', 'tSNE_2'])
        tsne_df['group'] = group_label
        return tsne_df

    group1_tsne = get_tsne_with_label(tsne_results[:len(group1)], "Group1")
    group2_tsne = get_tsne_with_label(tsne_results[len(group1):len(group1) + len(group2)], "Group2")
    group3_tsne = get_tsne_with_label(tsne_results[len(group1) + len(group2):len(group1) + len(group2) + len(group3)], "Group3")
    group4_tsne = get_tsne_with_label(tsne_results[len(group1) + len(group2) + len(group3):len(group1) + len(group2) + len(group3) + len(group4)], "Group4")

    combined_tsne_data = pd.concat([group1_tsne, group2_tsne, group3_tsne, group4_tsne])

    combined_tsne_data.to_csv('combined_tsne_results.csv', index=False)
