from mitofuncs.mitoevo import *
from mitofuncs.mitosignal import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

all_species_names = list_all_species_names_from_file_path()
all_species_genomes = {species: parse_mitochondrial_genome(f"data/species_files/{species}/{species}_mitochondrion.gb") for species in all_species_names}
all_species_gene_dataframes = {species: pd.read_csv(f"data/species_files/{species}/{species}_cleaned_gene_data.tsv", delimiter="\t").squeeze() for species in all_species_names}
all_species_gene_orders = {species: all_species_genomes[species]['Gene'].tolist() for species in all_species_names}

metric_name="damerau_levenshtein"
gene_orders_distance_matrix = []
for species_1 in all_species_names:
    row = []
    for species_2 in all_species_names:
        row.append(gene_order_distance(all_species_gene_orders[species_1], 
                                       all_species_gene_orders[species_2], 
                                       metric_name=metric_name))
    gene_orders_distance_matrix.append(row)
gene_orders_distance_dataframe = pd.DataFrame(gene_orders_distance_matrix, 
                                                 index=all_species_names, 
                                                 columns=all_species_names)
fig, ax = plt.subplots(figsize=(25,25))
hm = sn.heatmap(gene_orders_distance_dataframe,
                cmap='coolwarm', cbar_kws={"shrink": 1.0},
                square=True, 
                #vmin=0, vmax=1,
                #annot=True, fmt=".2g"
                )
fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
plt.title(f"Gene order distances calculated using {metric_name} metric")
plt.xlabel("species")
plt.ylabel("species")
plt.tight_layout()

gene_order_distances_folder_path = Path("results/gene_order_distances")
if not gene_order_distances_folder_path.exists():
    gene_order_distances_folder_path.mkdir(parents=True)

plt.savefig(f"results/gene_order_distances/gene_order_distances_heatmap.png", dpi=400)
plt.close()

MEGA_distance_df_string = convert_distance_dataframe_to_MEGA11_file_format(gene_orders_distance_dataframe, 
                                                                f"All gene order distances of the Current Dataset")
write_string_to_file(string=MEGA_distance_df_string, 
                        filepath=f"results/gene_order_distances/all_gene_order_distances.MEG")
gene_orders_distance_dataframe.to_csv(f"results/gene_order_distances/all_gene_order_distances.tsv", sep="\t")

print(gene_order_distance(all_species_gene_orders["Homo_sapiens"], 
                        all_species_gene_orders["Danio_rerio"], 
                        metric_name=metric_name))
print(gene_order_distance(all_species_gene_orders["Homo_sapiens"], 
                        all_species_gene_orders["Drosophila_melanogaster"],
                        metric_name=metric_name))