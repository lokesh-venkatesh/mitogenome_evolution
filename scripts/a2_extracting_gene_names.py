from mitofuncs.mito import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

all_species_names = list_all_species_names_from_file_path()

species_gene_data_dict = {species_name: parse_mitochondrial_genome(genbank_filepath=f"data/genbank_files/{species_name}/{species_name}_mitochondrion.gb")
                          for species_name in all_species_names}

def product_to_gene_series(species_gene_data_dict):
    all_data = [] # Iterate through the dictionary
    for species, df in species_gene_data_dict.items():
        if 'Product' in df.columns and 'Gene' in df.columns: # Ensure the required columns exist
            all_data.append(df[['Product', 'Gene']]) # Extract the Product-Gene mapping
        else:
            print(f"Warning: Missing 'Product' or 'Gene' columns in data for species: {species}")
    combined_df = pd.concat(all_data, ignore_index=True) # Concatenate all the dataframes and create the Series
    product_gene_series = combined_df.set_index('Product')['Gene']
    return product_gene_series

result_series = product_to_gene_series(species_gene_data_dict)
#result_series.to_csv("data/cleaned_gene_data/0_all_ACTUAL_raw_product_gene_mappings.tsv", sep='\t', header=True)

def remove_repeated_entries(series):
    unique_df = series.reset_index().drop_duplicates().set_index('Product') # Convert the Series to a DataFrame to drop duplicates based on both index and values
    unique_series = unique_df['Gene'] # Convert back to a Series
    return unique_series

final_result_series = remove_repeated_entries(result_series).sort_index()
final_result_series.to_csv("data/all_unique_product_gene_mappings.tsv", sep='\t', header=True)
final_result_series.to_csv("data/product_to_gene_mapping.tsv", sep='\t', header=True)