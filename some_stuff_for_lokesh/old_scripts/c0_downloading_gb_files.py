from pathlib import Path
from mitofuncs.mito import *
import pandas as pd
from Bio import Entrez

lokesh_email_ID = "lokesh.venkatesh@students.iiserpune.ac.in"
dataset_filepath = Path("data/current_species_dataset.tsv")


"""
def read_animal_ids(filepath):
    with open(filepath, 'r') as file:
        return [line.strip() for line in file.readlines()]

all_animal_mitogenomes_NCBI_IDs = read_animal_ids("data/animal_mitogenome_IDs.txt")

def fetch_organism_name(genbank_id, email):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
    record = handle.read()
    handle.close()
    for line in record.splitlines():
        if line.startswith("  ORGANISM"):
            return line.split()[1] + "_" + line.split()[2]
    return None

def create_animal_series(animal_ids, email):
    animal_dict = {}
    i = 0
    for genbank_id in animal_ids:
        print(i, genbank_id)
        organism_name = fetch_organism_name(genbank_id, email)
        if organism_name:
            animal_dict[genbank_id] = organism_name
            print(i,genbank_id, organism_name)
        i += 1
    return pd.Series(animal_dict)

animal_series = create_animal_series(all_animal_mitogenomes_NCBI_IDs, lokesh_email_ID)
animal_series.to_csv(dataset_filepath, sep="\t")
print("Dataset preparation is DONE")

animal_species_series = pd.read_csv(dataset_filepath, delimiter="\t", index_col=0).squeeze()
animal_species_series = pd.Series(animal_species_series.index, index=animal_species_series.values)
animal_species_series.to_csv(dataset_filepath, sep="\t")
"""


animal_species_series = pd.read_csv(dataset_filepath, delimiter="\t", index_col=0).squeeze()

genbank_files_folder_path = Path("genbank_files")
if not genbank_files_folder_path.exists(): 
    genbank_files_folder_path.mkdir(parents=True)

i=0
for animal, ID in animal_species_series.items():
    animal_folder_path = Path(f"genbank_files/{animal}")
    if not animal_folder_path.exists():
        animal_folder_path.mkdir(parents=True)

    gb_filepath = f"genbank_files/{animal}/{animal}_mitochondrion.gb"
    i+=1 
    print(i, animal)
    if not Path(gb_filepath).exists():
        fetch_and_save_genbank(genbank_id=ID, email=lokesh_email_ID, output_path=gb_filepath)
    else:
        print(f"Genbank file for {animal} already exists.")