from pathlib import Path
from mitofuncs.mito import *
import pandas as pd
from Bio import Entrez

lokesh_email_ID = "lokesh.venkatesh@students.iiserpune.ac.in"
dataset_filepath = Path("data/current_species_dataset.tsv")

animal_species_series = pd.read_csv(dataset_filepath, delimiter="\t", index_col=0).squeeze()

genbank_files_folder_path = Path("data/genbank_files")
if not genbank_files_folder_path.exists(): 
    genbank_files_folder_path.mkdir(parents=True)

i=0
for animal, ID in animal_species_series.items():
    animal_folder_path = Path(f"data/genbank_files/{animal}")
    if not animal_folder_path.exists():
        animal_folder_path.mkdir(parents=True)

    gb_filepath = f"data/genbank_files/{animal}/{animal}_mitochondrion.gb"
    i+=1 
    print(i, animal, ID)
    if not Path(gb_filepath).exists():
        fetch_and_save_genbank(genbank_id=ID, email=lokesh_email_ID, output_path=gb_filepath)
    else:
        print(f"Genbank file for {animal} already exists.")