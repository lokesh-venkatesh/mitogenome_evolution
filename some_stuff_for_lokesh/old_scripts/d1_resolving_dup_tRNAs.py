from mitofuncs.mito import *
import pandas as pd
from pathlib import Path

species_name = "Homo_sapiens"
gene_details_filepath = Path(f"genbank_files/{species_name}/{species_name}_cleaned_gene_data.tsv")

human_mitogenome = extract_genome_sequence(
                    Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb"))
Homo_sapiens_gene_df = pd.read_csv(gene_details_filepath, delimiter="\t").squeeze()
#print(Homo_sapiens_gene_df)

filtered_df = Homo_sapiens_gene_df.loc[Homo_sapiens_gene_df['Gene'].str.startswith('Trn')].set_index('Gene')

all_human_tRNA_genes = {}

for gene, row in filtered_df.iterrows():
    start = row['Start']-1
    stop = row['Stop']-1
    sequence = human_mitogenome[start:stop]
    all_human_tRNA_genes[gene] = sequence
'''
import subprocess
def get_rna_structure(sequence):
    result = subprocess.run(["RNAfold"], input=sequence, capture_output=True, text=True)
    return result.stdout

for key, value in all_human_tRNA_genes.items():
    tRNA_sequence = value.replace('T', 'U')
    print(key, get_rna_structure(tRNA_sequence))
'''

from Bio import SeqIO
def extract_specific_trna(genbank_file, trna_product_name):
    trna_entries = []  # List to store extracted tRNA data
    with open(genbank_file, "r") as gb_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type == "tRNA":
                    product = feature.qualifiers.get("product", [""])[0]
                    if product == trna_product_name:  # Match exact product name
                        codon_recognized = feature.qualifiers.get("codon_recognized", ["Not specified"])[0]
                        location = feature.location
                        trna_entries.append({
                            "Product": product,
                            "Codon Recognized": codon_recognized,
                            "Location": str(location)})
    return trna_entries

all_species = list_all_species_names_from_file_path()
for species_name in all_species:
    genbank_file = f"genbank_files/{species_name}/{species_name}_mitochondrion.gb"  # Replace with your actual file
    trna_product_name = "tRNA-Ser"  # Change to "tRNA-Ser" or any other tRNA
    trna_details = extract_specific_trna(genbank_file, trna_product_name)
    if trna_details:
        print(f"\nDetails for {trna_product_name}: where species is {species_name}")
        for i, entry in enumerate(trna_details, 1):
            print(f"  {i}. Location: {entry['Location']}, Codon Recognized: {entry['Codon Recognized']}")
    else:
        print(f"No entries found for {trna_product_name}.")