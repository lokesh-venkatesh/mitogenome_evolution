from mitofuncs.mitoevo import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

all_species_names = list_all_species_names_from_file_path()
all_cleaned_gene_dataframes = {species: pd.read_csv(Path(f"genbank_files/{species}/{species}_cleaned_gene_data.tsv"), 
                                                    delimiter="\t") for species in all_species_names}
all_gene_orders = {species: list(all_cleaned_gene_dataframes[species]['Gene']) for species in all_species_names}

def combine_gene_orders_to_dataframe(gene_orders_dict):
    combined_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in gene_orders_dict.items()]))
    return combined_df

folder_path = Path("results/visualising_translocations")
if not folder_path.exists():
    folder_path.mkdir(parents=True)

combined_gene_orders_df = combine_gene_orders_to_dataframe(all_gene_orders)
#combined_gene_orders_df = combined_gene_orders_df[['Arabidopsis_thaliana','Saccharomyces_cerevisiae','Caenorhabditis_elegans','Drosophila_melanogaster','Homo_sapiens','Mus_musculus','Xenopus_tropicalis','Danio_rerio']]
combined_gene_orders_df.to_csv("results/visualising_translocations/0_combined_gene_orders.tsv", sep="\t", index=False)


def gene_mapping(gene_order_list_1, gene_order_list_2):
    mapping = []
    for gene in gene_order_list_1:
        indices = [i for i, x in enumerate(gene_order_list_2) if x == gene]  # Find all indices of the gene in list2
        mapping.append(indices if indices else [-1])  # Append -1 if the gene is not found
    return mapping

LMAO = 1
for species_1, gene_order_list_1 in all_gene_orders.items():
    print(f"{species_1}, aka {LMAO}th genome going on right now, out of 61")
    LMAO += 1
    for species_2, gene_order_list_2 in all_gene_orders.items():
        file_path = folder_path / Path(f"{species_1}_and_{species_2}_gene_orders_visualised.png")
        if not os.path.exists(file_path):
            mapping = gene_mapping(gene_order_list_1, gene_order_list_2) # Get the mapping of genes from gene_order_list_1 to gene_order_list_2

            x_vals = []
            y_vals = []

            for i, indices in enumerate(mapping):
                for index in indices:
                    if index != -1:  # Ignore unmapped genes
                        x_vals.append(i)
                        y_vals.append(index)

            plt.figure(figsize=(12, 12)) # Plotting the linear mapping
            plt.scatter(x_vals, y_vals, color="b")
            plt.plot(x_vals, y_vals, linestyle="--", alpha=0.5, color="gray") # Plot the linear mapping
            plt.xticks(range(len(gene_order_list_1)), gene_order_list_1, rotation=90) # Set axis labels and title
            plt.yticks(range(len(gene_order_list_2)), gene_order_list_2, rotation=0)
            plt.xlabel(f"Order of genes in {species_1}")
            plt.ylabel(f"Order of genes in {species_2}")
            plt.title(f"Gene Orders of {species_1} versus {species_2}")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(folder_path / Path(f"{species_1}_and_{species_2}_gene_orders_visualised.png"),dpi=300)
            plt.close()