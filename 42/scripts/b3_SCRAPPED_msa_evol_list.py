from mitofuncs.mitoevo import *
from mitofuncs.mitosignal import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
import os
import subprocess
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import AlignIO

all_species_names = list_all_species_names_from_file_path()

all_genomes_as_one_fasta_string = return_all_mitogenome_sequences_as_one_fasta(all_species_names)
with open("data/all_genomes.fasta", "w") as f:
    f.write(all_genomes_as_one_fasta_string)

'''
import subprocess
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import os

# CONFIG: Define tools and parameters
MUSCLE_EXEC = "C:/Users/lokes/Desktop/muscle-win64.v5.3.exe"  # Ensure this is in your system PATH
TREE_METHOD = "nj"  # Use "nj" (Neighbor-Joining) or "upgma"

# INPUT FILE (Modify as needed)
INPUT_FASTA = "data/all_genomes.fasta"
ALIGNED_FASTA = "data/all_genomes_aligned.fasta"
TREE_FILE = "data/all_genomes_MSA_tree.nwk"

# Step 1: Run MSA using MUSCLE
def run_muscle(input_fasta, output_fasta):
    """Runs MUSCLE for multiple sequence alignment."""
    subprocess.run([MUSCLE_EXEC, "-align", input_fasta, "-output", output_fasta], check=True)

# Step 2: Compute Distance Matrix from Aligned Sequences
def compute_distance_matrix(aligned_fasta):
    """Computes pairwise sequence distances to build a distance matrix."""
    alignment = AlignIO.read(aligned_fasta, "fasta")
    num_sequences = len(alignment)

    # Compute pairwise distances based on sequence identity
    matrix = [[0] * num_sequences for _ in range(num_sequences)]
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            matches = sum(1 for a, b in zip(alignment[i].seq, alignment[j].seq) if a == b)
            distance = 1 - (matches / len(alignment[i].seq))  # Normalized distance
            matrix[i][j] = matrix[j][i] = distance

    names = [record.id for record in alignment]
    return DistanceMatrix(names, matrix)

# Step 3: Build Phylogenetic Tree
def build_tree(distance_matrix, method="nj"):
    """Builds a phylogenetic tree using Neighbor-Joining (NJ) or UPGMA."""
    constructor = DistanceTreeConstructor()
    if method == "upgma":
        tree = constructor.upgma(distance_matrix)
    else:
        tree = constructor.nj(distance_matrix)
    return tree

# Step 4: Extract Evolutionary Order from the Tree
def extract_evolutionary_order(tree):
    """Extracts an ordered list of sequences from the phylogenetic tree."""
    return [leaf.name for leaf in tree.get_terminals()]

# Step 5: Full Pipeline Execution
def process_fasta(input_fasta):
    print(f"Running MSA on {input_fasta}...")
    run_muscle(input_fasta, ALIGNED_FASTA)

    print("Computing distance matrix...")
    distance_matrix = compute_distance_matrix(ALIGNED_FASTA)

    print(f"Building {TREE_METHOD} tree...")
    tree = build_tree(distance_matrix, TREE_METHOD)

    print("Saving tree to file...")
    Phylo.write(tree, TREE_FILE, "newick")

    print("Extracting evolutionary order...")
    evolutionary_order = extract_evolutionary_order(tree)

    return evolutionary_order

if __name__ == "__main__":
    evolutionary_order = process_fasta(INPUT_FASTA)
    print("\nEvolutionary Order of Sequences:")
    print("\n".join(evolutionary_order))'
'''