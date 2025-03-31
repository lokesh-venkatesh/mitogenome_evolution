import mitofuncs.mitoevo as mito
from pathlib import Path
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import time

lokesh_email_ID = "lokesh.venkatesh@students.iiserpune.ac.in"
dataset_filepath = Path("data/animal_mitochondrial_genomes.csv")

# Assuming the CSV file has two columns: 'animal' and 'ID'
animal_species_df = pd.read_csv(dataset_filepath, header=0)
animal_species_df['Genome Name'] = animal_species_df['Genome Name'].str.extract(r'^(.*?) mitochondrion, complete genome$')[0]
print(animal_species_df.head())

species_files_folder_path = Path("data/species_files")
if not species_files_folder_path.exists(): 
    species_files_folder_path.mkdir(parents=True)

start_time = time.time()

for index, row in animal_species_df.iterrows():
    animal = row['Genome Name']
    ID = row['NCBI Accession']

    animal = animal.replace(" ", "_")
    animal_folder_path = Path(f"data/species_files/{animal}")
    if not animal_folder_path.exists():
        animal_folder_path.mkdir(parents=True)
    gb_filepath = f"data/species_files/{animal}/{animal}_mitochondrion.gb"
    print(animal, ID)
    if not Path(gb_filepath).exists():
        mito.fetch_and_save_genbank(genbank_id=ID, email=lokesh_email_ID, output_path=gb_filepath)
    else:
        print(f"Genbank file for {animal} already exists.")

end_time = time.time()
print(f"Time taken to run the script: {end_time - start_time} seconds")