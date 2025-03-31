import argparse
import os

# Argument parser for dynamic output directory
parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", type=str, default="data")  # Default to "data" if not provided
args = parser.parse_args()

# Ensure directory exists
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# Define output file path
output_file = os.path.join(output_dir, "animal_mitochondrial_genomes.csv")

print(f"Saving results to: {output_file}")


import csv
from Bio import Entrez
import time

# Set your email (NCBI requires this)
Entrez.email = "lokesh.venkatesh@students.iiserpune.ac.in"

# Define query for complete animal mitochondrial reference genomes
query = '("mitochondrion"[All Fields] AND "complete genome"[Title]) AND "RefSeq"[filter] AND "Animals"[Organism]'

# Step 1: Get total number of records dynamically
print("Querying NCBI for total count of animal mitochondrial genomes...")
handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)  # Get only count
record = Entrez.read(handle)
handle.close()
total_records = int(record["Count"])  # Get total available records

print(f"Total records available: {total_records}")

# Step 2: Fetch all IDs in batches
batch_size = 5000  # NCBI's safe batch limit
ncbi_ids = []

print("Fetching genome IDs...")
for start in range(0, total_records, batch_size):
    print(f"Fetching records {start} to {min(start + batch_size, total_records)}...")
    handle = Entrez.esearch(db="nucleotide", term=query, retstart=start, retmax=batch_size)
    record = Entrez.read(handle)
    handle.close()
    ncbi_ids.extend(record["IdList"])
    time.sleep(0.5)  # Prevent overwhelming NCBI servers

print(f"Total genomes fetched: {len(ncbi_ids)}")

# Step 3: Fetch metadata in batches
data = []
batch_size_summary = 500  # NCBI has limits on esummary batch sizes
print("Fetching metadata for genomes...")

for start in range(0, len(ncbi_ids), batch_size_summary):
    print(f"Processing metadata for records {start} to {min(start + batch_size_summary, len(ncbi_ids))}...")
    batch_ids = ",".join(ncbi_ids[start:start + batch_size_summary])
    
    handle = Entrez.esummary(db="nucleotide", id=batch_ids)
    summary_records = Entrez.read(handle)
    handle.close()
    
    for summary in summary_records:
        accession = summary["AccessionVersion"]
        title = summary["Title"]
        data.append([accession, title])
        if len(data) % 100 == 0:
            print(f"Processed {len(data)} genomes")

    time.sleep(0.5)  # Prevent hitting rate limits

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["NCBI Accession", "Genome Name"])
    writer.writerows(data)

print(f"Saved results to {output_file}")

#scp C:/Users/lokes/Desktop/mitoevo/scripts/a0_obtaining_dataset.py madhu@192.168.1.152:/storage/madhu/lokesh/42_jobs/
#scp C:/Users/lokes/Desktop/mitoevo/42_jobs/anim_genomes_dwnld.sh madhu@192.168.1.152:/storage/madhu/lokesh/42_jobs/