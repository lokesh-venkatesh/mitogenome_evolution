import mitofuncs.mitoevo as mito
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

all_species_names = mito.list_all_species_names_from_file_path()

species_gene_data_dict = {species_name: mito.parse_mitochondrial_genome(genbank_filepath=f"data/species_files/{species_name}/{species_name}_mitochondrion.gb")
                          for species_name in all_species_names}

def product_to_gene_series(species_gene_data_dict):
    """Converts all of the items in the product column of the gene data into a pandas series
    with products as index column and gene abbreviations as the value column"""
    all_data = []
    for species, df in species_gene_data_dict.items():
        if 'Product' in df.columns and 'Gene' in df.columns: 
            all_data.append(df[['Product', 'Gene']]) 
        else:
            print(f"Warning: Missing 'Product' or 'Gene' columns in data for species: {species}")
    combined_df = pd.concat(all_data, ignore_index=True)
    product_gene_series = combined_df.set_index('Product')['Gene']
    return product_gene_series

result_series = product_to_gene_series(species_gene_data_dict)

def remove_repeated_entries(series):
    unique_df = series.reset_index().drop_duplicates().set_index('Product') 
    unique_series = unique_df['Gene']
    return unique_series

final_result_series = remove_repeated_entries(result_series).sort_index()
final_result_series.to_csv("data/all_unique_products_to_genes.tsv", sep='\t', header=True)

"""NOTE: THE BELOW .TSV FILE MUST BE MANUALLY EDITED FOR FURTHER WORK SO THAT THE GENE DATA IS CLEANED UP
THE WAY YOU WANT IT TO BE. THIS IS BECAUSE THE GENE DATA IS NOT STANDARDIZED ACROSS ALL SPECIES AND
THEREFORE NEEDS TO BE CLEANED UP MANUALLY."""

final_result_series.to_csv("data/shortlisted_products_to_genes.tsv", sep='\t', header=True)