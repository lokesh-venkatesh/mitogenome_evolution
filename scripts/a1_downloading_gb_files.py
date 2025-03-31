from pathlib import Path
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
import time

def fetch_and_save_genbank(genbank_id, email, output_path):
    """
    Fetches a GenBank file for a given ID and saves it to a specified location.
    
    Parameters:
    genbank_id (str): The GenBank accession number or GI number.
    email (str): Email ID required for NCBI Entrez API usage.
    output_path (str): Full path including filename to save the GenBank file.
    """
    Entrez.email = email  # Set the email ID for NCBI access
    try:
        # Fetch the GenBank record
        with Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gbwithseq", retmode="text") as handle:
            records = list(SeqIO.parse(handle, "genbank"))  # Read all records before closing handle
        # Save the file
        if records:
            with open(output_path, "w") as output_file:
                SeqIO.write(records, output_file, "genbank")  # Write in GenBank format
            print(f"GenBank file saved successfully at: {output_path}")
        else: print(f"ERROR: No records found for {genbank_id}")
    except Exception as e:
        print(f"ERROR FETCHING THE GENBANK RECORD: {e}")

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
        fetch_and_save_genbank(genbank_id=ID, email=lokesh_email_ID, output_path=gb_filepath)
    else:
        print(f"Genbank file for {animal} already exists.")

end_time = time.time()
print(f"Time taken to run the script: {end_time - start_time} seconds")