from mitofuncs.mitoevo import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

folder_path = Path("results/gene_orders")
second_folder_path = Path("genbank_files")

if not folder_path.exists():
    folder_path.mkdir(parents=True)
    print(f"Folder '{folder_path}' created.")
else: 
    print(f"Folder '{folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()

third_path = Path("results/cleaned_gene_data")

if not third_path.exists():
    third_path.mkdir(parents=True)
    print(f"Folder '{third_path}' created.")
else: 
    print(f"Folder '{third_path}' already exists.")


species_gene_data_dict = {species_name: parse_mitochondrial_genome(genbank_filepath=f"genbank_files/{species_name}/{species_name}_mitochondrion.gb") 
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
result_series.to_csv("results/cleaned_gene_data/0_all_ACTUAL_raw_product_gene_mappings.tsv", sep='\t', header=True)

def remove_repeated_entries(series):
    unique_df = series.reset_index().drop_duplicates().set_index('Product') # Convert the Series to a DataFrame to drop duplicates based on both index and values
    unique_series = unique_df['Gene'] # Convert back to a Series
    return unique_series

final_result_series = remove_repeated_entries(result_series).sort_index()
final_result_series.to_csv("results/cleaned_gene_data/0_all_unique_product_gene_mappings.tsv", sep='\t', header=True)


repeated_keys = final_result_series.index[final_result_series.index.duplicated()].unique()
filtered_data = final_result_series[~final_result_series.index.isin(repeated_keys)]
filtered_data.to_csv("results/cleaned_gene_data/0_all_FINAL_unique_product_gene_mappings_filtered.tsv", sep='\t', header=True)
final_result_series.to_csv("results/cleaned_gene_data/0_all_FINAL_unique_product_gene_mappings_filtered.tsv", sep='\t', header=True)


"""
WE USE THE CREATED '0_all_FINAL_unique_product_gene_mappings_filtered.tsv' FILE,
TO CREATE THE 'product_to_gene_mapping' DICTIONARY

FOR THIS, THIS CODE FILE WAS RUN TWICE.
"""

def update_gene_column(gene_data_df, product_to_gene_mapping):
    if 'Product' not in gene_data_df.columns or 'Gene' not in gene_data_df.columns:
        raise ValueError("Input DataFrame must contain 'Product' and 'Gene' columns.") # Update the 'Gene' column based on the mapping Series
    gene_data_df['Gene'] = gene_data_df['Product'].map(product_to_gene_mapping).fillna(gene_data_df['Gene'])
    return gene_data_df

product_to_gene_mapping = pd.read_csv("results/cleaned_gene_data/0_product_to_gene_mapping.tsv", sep='\t', header=0).squeeze()
product_to_gene_mapping.columns = ['Product', 'Gene']
mapping_series = pd.Series(product_to_gene_mapping['Gene'].values, index=product_to_gene_mapping['Product'])
mapping_series = product_to_gene_mapping

for species_name in all_species_names:
    gb_filepath = f"genbank_files/{species_name}/{species_name}_mitochondrion.gb"
    gene_data = parse_mitochondrial_genome(genbank_filepath=gb_filepath)
    gene_data = update_gene_column(gene_data, mapping_series)
    gene_data.to_csv(folder_path / f"{species_name}_gene_data.tsv", sep="\t", index=False)
    gene_data.to_csv(second_folder_path / f"{species_name}" / f"{species_name}_cleaned_gene_data.tsv", sep="\t", index=False)
    gene_data.to_csv("results/cleaned_gene_data/"+f"{species_name}_cleaned_gene_data.tsv", sep="\t", index=False)
