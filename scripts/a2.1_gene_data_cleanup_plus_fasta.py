"""
THIS FILE IS MEANT TO USE THE produce_to_gene_mapping.tsv FILE AND CLEAN UP ALL OF THE GENE DATA,
WHICH IS THEN SAVED IN THE CORRESPONDING FOLDERS OF THE GENBANK FILES FOR EACH SPECIES.
"""


from mitofuncs.mito import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

"""
WE USE THE CREATED '0_all_FINAL_unique_product_gene_mappings_filtered.tsv' FILE,
TO CREATE THE 'product_to_gene_mapping' DICTIONARY IN THE 'product_to_gene_mapping.tsv' FILE
"""

all_species_names = list_all_species_names_from_file_path()

def update_gene_column(gene_data_df, product_to_gene_mapping):
    if 'Product' not in gene_data_df.columns or 'Gene' not in gene_data_df.columns:
        raise ValueError("Input DataFrame must contain 'Product' and 'Gene' columns.") # Update the 'Gene' column based on the mapping Series
    gene_data_df['Gene'] = gene_data_df['Product'].map(product_to_gene_mapping).fillna(gene_data_df['Gene'])
    return gene_data_df

product_to_gene_mapping = pd.read_csv("data/product_to_gene_mapping.tsv", sep='\t', header=0).squeeze()
mapping_series = product_to_gene_mapping.set_index("Product")["Gene"].to_dict()

for species_name in all_species_names:
    gb_filepath = f"data/genbank_files/{species_name}/{species_name}_mitochondrion.gb"
    gene_data = parse_mitochondrial_genome(genbank_filepath=gb_filepath)
    gene_data = update_gene_column(gene_data, mapping_series)
    gene_data.to_csv(Path("data/genbank_files") / f"{species_name}" / f"{species_name}_cleaned_gene_data.tsv", sep="\t", index=False)
    genome_fasta_string = f">{species_name}:\n{extract_genome_sequence(genbank_filepath=gb_filepath)}"
    write_string_to_file(string=genome_fasta_string, filepath=f"data/genbank_files/{species_name}/{species_name}_mitogenome.fasta")