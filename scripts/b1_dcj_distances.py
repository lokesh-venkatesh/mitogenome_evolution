from pathlib import Path
from mitofuncs.mito import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

all_species_names = list_all_species_names_from_file_path()
all_species_gene_data = {species: pd.read_csv(f"data/genbank_files/{species}/{species}_cleaned_gene_data.tsv", 
                                              delimiter="\t").squeeze() for species in all_species_names}
all_species_gene_order_lists = {species: all_species_gene_data[species]['Gene'].tolist() 
                                for species in all_species_names}

all_genomes_DCJ_breakpoint_distance_matrix = []
for species_1 in all_species_names:
    gene_order_1 = all_species_gene_order_lists[species_1]
    row = []
    for species_2 in all_species_names:
        gene_order_2 = all_species_gene_order_lists[species_2]
        row.append(double_cut_and_join_model_breakpoint_distance(gene_order_1, gene_order_2))
    all_genomes_DCJ_breakpoint_distance_matrix.append(row)
    print(species_1)

all_genomes_DCJ_breakpoint_distance_df = pd.DataFrame(all_genomes_DCJ_breakpoint_distance_matrix, 
                                                   index=all_species_names, columns=all_species_names)
plt.figure(figsize=(12,12))
hm = sn.heatmap(all_genomes_DCJ_breakpoint_distance_df,
                cmap='coolwarm',
                square=True, 
                vmin=0)
plt.title(f"Heatmap of all DCJ_breakpoint distances for all species' gene orders taken pairwise")
plt.xlabel("Species")
plt.ylabel("Species")
plt.tight_layout()

DCJ_breakpoint_folder_path = Path("results/DCJ_breakpoint")
if not DCJ_breakpoint_folder_path.exists():
    DCJ_breakpoint_folder_path.mkdir(parents=True)

plt.savefig(f"results/DCJ_breakpoint/all_genomes_DCJ_breakpoint_distance_df.png", dpi=400)
all_genomes_DCJ_breakpoint_distance_df.to_csv(f"results/DCJ_breakpoint/all_genomes_DCJ_breakpoint_distance_df", sep="\t")

MEGA11_DCJ_breakpoint_Distances_String = convert_dataframe_to_MEGA11_file_format(all_genomes_DCJ_breakpoint_distance_df, "DCJ_breakpoint Distances for all Gene Orders")
write_string_to_file(MEGA11_DCJ_breakpoint_Distances_String, "results/DCJ_breakpoint/all_genomes_DCJ_breakpoint_distances.MEG")