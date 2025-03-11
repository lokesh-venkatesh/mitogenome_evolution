import os
import Bio
import random
import itertools
import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import Entrez,SeqIO
import time
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree
import ete3
from ete3 import Tree

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

def extract_genome_sequence(genbank_filepath):
    """READS THE GENBANK FILEPATH AND RETURNS THE GENOME SEQUENCE CONTAINED IN THE .gb FILE"""
    root_dir = Path(__file__).parent.parent
    file_path = root_dir / genbank_filepath # Ensure the file exists before proceeding
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    from Bio import SeqIO # Proceed with loading the file
    record = SeqIO.read(file_path, "genbank")
    return str(record.seq)

def parse_mitochondrial_genome(genbank_filepath):
    """Parse a mitochondrial genome GenBank (.gb) file to extract information about CDS and RNA sequences.
    Parameters: gb_file (str): Path to the GenBank file.
    Returns: pd.DataFrame: A DataFrame containing gene name, gene symbol, start position, stop position, strand, sequence, and product."""
    record = SeqIO.read(genbank_filepath, "genbank")
    genome_seq = str(record.seq)  # Get the full genome sequence as a string
    genome_length = len(genome_seq)  # Total genome length
    data = []
    for feature in record.features:
        if feature.type in ["CDS", "rRNA", "tRNA"]:  
            gene_symbol = feature.qualifiers.get("locus_tag", [None])[0]
            product = feature.qualifiers.get("product", [None])[0]
            gene = feature.qualifiers.get("gene", [None])[0]

            start = int(feature.location.start) + 1
            stop = int(feature.location.end)

            strand = "+" if feature.location.strand == 1 else "-"

            if start > stop:
                sequence = genome_seq[start-1:] + genome_seq[:stop]
            else:
                sequence = genome_seq[start-1:stop]

            data.append([gene, product, start, stop, strand])
    df = pd.DataFrame(data, columns=["Gene", "Product", "Start", "Stop", "Strand"
                                     ])
    """ MAYBE CONSIDER ADDING 'sequence' AS WELL TO THE ABOVE LIST NAMED 'data' """
    return df

def list_all_species_names_from_file_path(file_path="data/current_species_dataset.tsv"):
    """RETURNS THE SPECIES NAMES IN THE GENBANK FOLDERS' DIRECTORY"""
    try: # Get all files in the directory
        with open(file_path, 'r') as file: # If the file exists, load the output
            animal_species_series = pd.read_csv(file_path, sep="\t", index_col=0).squeeze()
            all_species_names = list(animal_species_series.keys())
        return all_species_names
    except FileNotFoundError:
        print(f"Error: The directory '{file_path}' does not exist.")
        return []
    except PermissionError:
        print(f"Error: Permission denied to access '{file_path}'.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []
    
def return_all_mitochondrion_genome_sequences_as_dict(species_names, folder_path="genbank_files"):
    """RETURNS A DICTIONARY CONTAINING THE GENOME SEQUENCES OF ALL SPECIES IN THE SPECIES_NAMES LIST"""
    genome_sequences = {}
    for species_name in species_names:
        root_dir = Path(__file__).parent.parent
        genbank_filepath = os.path.join(root_dir, folder_path, species_name, f"{species_name}_mitochondrion.gb")
        genome_sequences[species_name] = extract_genome_sequence(genbank_filepath)
    return genome_sequences

def return_genome_nucl_frequencies_as_dict(genome_sequence):
    """RETURNS A DICTIONARY CONTAINING THE NUCLEOTIDE FREQUENCIES OF THE INPUT-GENOME SEQUENCE"""
    try:
        frequencies = {"A": 0, "T": 0, "C": 0, "G": 0}
        for nucleotide in genome_sequence.upper():  # Convert to uppercase to handle case-insensitivity
            if nucleotide in frequencies:  # Check if it's a valid nucleotide
                frequencies[nucleotide] += 1
        
        for nucleotide in frequencies.keys():
            frequencies[nucleotide] = frequencies[nucleotide]/len(genome_sequence)
        return frequencies
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def sliding_window_nucleotide_frequencies(sequence, window_size=100, step_size=1):
    """DOES A MOVING WINDOW ANALYSIS OF THE NUCLEOTIDE FREQUENCIES ALONG A GIVEN GENOME FOR A SPECIFIED WINDOW AND STEP SIZE
    RETURNS THE POSITION INDICES FOLLOWED BY THE NUCLEOTIDE_FREQUENCY LISTS AS A DICTIONARY"""
    nucleotides = ['A', 'T', 'G', 'C']
    positions = []
    frequencies = {nuc: [] for nuc in nucleotides}

    circular_sequence = sequence+sequence
    for i in range(0, len(sequence) + 1, step_size):
        window = circular_sequence[i:i + window_size]
        positions.append(i)
        total_bases = len(window)
        for nuc in nucleotides:
            frequencies[nuc].append(window.count(nuc) / total_bases)
    return positions, frequencies

def generate_kmer_motif_list(k, folder_path="data/motif_lists"):
    """GENERATES ALL POSSIBLE NUCLEOTIDES OF LENGTH==k AND RETURNS A PYTHON LIST, 
    OR JUST FETCHES IT FROM THE DUMP FOLDER"""

    def generate_nucleotide_list():
        """RETURNS A LIST OF THE FOUR NUCLEOTIDES"""
        return ["A","T","G","C"]
    
    def extend_iterative_motif_list(kmer_list):
        """TAKES EACH KMER IN THE GIVEN LIST, ADDS AN EXTRA NUCLEOTIDE TO ITS END FOR EACH SUCH NUCLEOTIDE, AND THEN RETURNS THE FINAL LIST OF ALL SUCH (K+1)MERS AS A PYTHON LIST"""
        extended_kmer_list = []
        nucleotides = generate_nucleotide_list()
        for kmer in kmer_list:
            for nucl in nucleotides:
                extended_kmer_list.append(kmer+nucl)
        return extended_kmer_list
    
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_path = os.path.join(folder_path, f"motif_list_{k}mers") # Define the path to the output file
    if os.path.exists(file_path): # Check if the output file exists
        with open(file_path, 'r') as file: # If the file exists, load the output
            output = [line.strip() for line in file]
        print(f"Loaded output for {k} from {file_path}.")
    else: # If the file does not exist, generate and save the output
        output = generate_nucleotide_list()
        for i in range(k-1):
            output = extend_iterative_motif_list(output)
        with open(file_path, 'w') as file:
            for item in output:
                file.write(f"{item}\n")
        print(f"Generated and saved output for {k}.")
    return output

def motif_observation_in_genome(genome_sequence,motif):
    """RETURNS THE NUMBER OF OBSERVED INSTANCES OF A MOTIF STRING IN THE GIVEN GENOME STRING"""
    motif_count = 0
    for i in range(len(genome_sequence)-len(motif)+1):
        if genome_sequence[i:i+len(motif)] == motif:
            motif_count+=1
    return motif_count

def motif_expectation_in_genome(genome_sequence,motif,nucl_freq):
    """RETURNS THE EXPECTED NUMBER OF OCCURENCES OF A MOTIF IN A SEQUENCE USING THE EXPECTATION FORMULA
    NOTE THAT NUCLEOTIDE FREQUENCIES MUST BE THAT OF THE WHOLE GENOME AND NOT JUST OF THAT SEQUENCE!!!"""
    expectation = 1
    for i in range(len(motif)):
        expectation = expectation*nucl_freq[motif[i]]
    return expectation*len(genome_sequence)

def gc_content(motif):
    """CALCULATES THE NUMBER OF GC NUCLEOTIDES IN THE GIVEN MOTIF"""
    return sum([1 for x in motif if x in ["G","C"]])

def sort_kmer_vector_by_GC_content(kmer_vector, species_name, k):
    """TAKES THE GIVEN KMER VECTOR AND RETURNS THE SAME LIST BUT SORTED BY INCREASING GC CONTENT
    KNOW THAT IT DOES NOT TAKE INTO ACCOUNT THE NUCLEOTIDE POSITIONS BETWEEN TWO MOTIFS WITH THE SAME GC CONTENT"""
    sorted_series = kmer_vector.sort_index(key=lambda x: x.map(gc_content))
    sorted_series.name = f"{species_name} {k}-mer OE values GC sorted"
    return sorted_series
    '''
    motif_list = list(kmer_vector.keys())
    sorted_motif_list = sorted(motif_list, key=gc_content)
    sorted_OE_ratios = []
    for motif in sorted_motif_list:
        sorted_OE_ratios.append(kmer_vector[motif])
    return pd.Series(sorted_OE_ratios,index=sorted_motif_list,name=f"{species_name} {k}-mer OE values GC sorted")
    '''

def generate_full_kmer_vector(k, species_name, genome_sequence, nucl_freq, file_name, folder_path="data/whole_genome_kmer_vectors"):
    """GENERATES THE FULL SIZE K-MER VECTOR FOR THE GIVEN GENOME SEQUENCE"""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_path = os.path.join(folder_path, file_name) # Define the path to the output file
    if os.path.exists(file_path): # Check if the output file exists
        output = pd.read_csv(file_path, delimiter="\t", index_col=0).squeeze()
        print(f"Loaded output for {k} from {file_path}.")
    else: # If the file does not exist, generate and save the output
        motif_list = generate_kmer_motif_list(k)
        OE_list = []
        for motif in motif_list:
            observed = motif_observation_in_genome(genome_sequence,motif)
            expected = motif_expectation_in_genome(genome_sequence,motif,nucl_freq)
            OE_ratio = observed/expected
            OE_list.append(OE_ratio) #Use this for dropping all but up to 3 decimal places: "%.3f"%
        kmer_vector = pd.Series(OE_list,index=motif_list,name=f"{species_name} {k}-mer full OE values")
        output = sort_kmer_vector_by_GC_content(kmer_vector, species_name, k)
        output.to_csv(file_path, sep="\t")
        print(f"Generated and saved output for {k}.")
    return output

def RC(motif):
    """RETURNS THE REVERSE COMPLEMENT OF A GIVEN SEQUENCE"""
    rc_dict = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    rc_motif = ""
    for i in range(len(motif)-1,-1,-1):
        if motif[i] in list(rc_dict.keys()):
            rc_motif += rc_dict[motif[i]]
        else:
            rc_motif += motif[i]
    return rc_motif

def generate_reduced_kmer_vector(raw_full_vector):
    """GENERATES THE NOTE REDUCED K-MER VECTOR FOR THE GIVEN GENOME SEQUENCE"""
    try: # If the file does not exist, generate and save the output
        full_vector = raw_full_vector.to_dict()
        all_kmers = list(full_vector.keys())
        exhausted_kmers = []
        reduced_kmer_dict = {}
        for kmer in all_kmers:
            if kmer in exhausted_kmers or RC(kmer) in exhausted_kmers: 
                continue
            else: 
                reduced_kmer_dict[kmer+"/"+RC(kmer)] = float((full_vector[kmer]+full_vector[RC(kmer)])/2)  # "%.3f"%
                exhausted_kmers.append(kmer)
        output = pd.Series(list(reduced_kmer_dict.values()),index=list(reduced_kmer_dict.keys()),name="reduced OE values")
        return output
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def pcc(series_1,series_2):
    """RETURNS THE PEARSON CORRELATION COEFFICIENT BETWEEN TWO MOTIF VECTORS"""
    if len(series_1)!=len(series_2):
        raise ValueError("Both series must have the same length")
    r = series_1.corr(series_2)
    return r

def pearson_correlation_coefficient(series_1,series_2):
    """RETURNS THE PEARSON CORRELATION COEFFICIENT BETWEEN TWO MOTIF VECTORS"""
    if len(series_1)!=len(series_2):
        raise ValueError("Both series must have the same length")
    r = series_1.corr(series_2)
    return r

def randomise(genome_seq):
    """SCRAMBLES THE GIVEN GENOME USING THE NUCLEOTIDE FREQUENCIES OBTAINED FROM THE ORIGINAL GENOME"""
    sequence_list = list(genome_seq)
    random.shuffle(sequence_list)
    return ''.join(sequence_list)

def generate_correlation_matrix(dictionary_of_vectors, r_type):
    """DOES ALL NC2 CORRELATIONS BETWEEN PAIRWISE-CHOSEN VECTORS FROM THE INPUT DICTIONARY
    CHOOSE EITHER r_type="pcc" pr r_type="modified_pcc" for this function!!! """
    species_names = list(dictionary_of_vectors.keys())
    matrix = []
    for vector_1 in dictionary_of_vectors.values():
        matrix_row = []
        for vector_2 in dictionary_of_vectors.values():
            if r_type=="pcc":
                r = pcc(vector_1, vector_2)
            elif r_type=="modified_pcc":
                r = pearson_correlation_coefficient(vector_1, vector_2)
            matrix_row.append(r)
        matrix.append(matrix_row)
    correlation_dataframe = pd.DataFrame(matrix, index=species_names, columns=species_names)
    return correlation_dataframe

def generate_distance_matrix(dictionary_of_vectors, r_type):
    """CALCULATES THE DISTANCE MATRIX FROM THE CORRELATION MATRIX GENERATED
    CHOOSE EITHER r_type="pcc" pr r_type="modified_pcc" for this function!!! """
    corrln_df = generate_correlation_matrix(dictionary_of_vectors, r_type)
    corrln_df = (corrln_df+1)/2
    distance_df = -np.log10(corrln_df)
    return distance_df

def convert_distance_dataframe_to_MEGA11_file_format(distance_df,title):
    """TAKES A DISTANCE MATRIX (STORED AS A DATAFRAME) AND CONVERTS IT INTO A STRING IN THE EXACT .MEG FILE FORMAT WITHOUT ANY ERRORS"""
    n_rows, n_cols = distance_df.shape
    file_string = f"#mega\n!Title {title} distance matrix;\n!Format DataType=Distance;\n\n"
    animal_species_list = distance_df.columns
    for animal in animal_species_list:
        file_string += f"#{animal}\n"
    file_string += "\n\n"
    row_ind_string = "\t[0"
    for i in range(1,n_rows):
        row_ind_string+=f"\t{i}"
    file_string+= f"{row_ind_string}]\n"

    for i in range(n_rows):
        new_row = f"[{i}]\t"
        for j in range(i): 
            new_row += f"{distance_df.iloc[i,j]}\t"
        file_string+=f"{new_row}\n"
    return file_string

def top_X_truncated_motifs_vector_pearsons_correlation(vector_1, vector_2, X):
    """TAKES ONLY THE TOP X VALUES OF EACH VECTOR,
    CONSTRUCTS THE UNION VECTOR FOR EACH, AND THEN RETURNS THE CORRELATION SCORE"""
    top_indices_1 = vector_1.nlargest(X).index
    top_indices_2 = vector_2.nlargest(X).index
    combined_indices = top_indices_1.union(top_indices_2)
    reduced_vector_1 = vector_1[combined_indices]
    reduced_vector_2 = vector_2[combined_indices]
    r = pcc(reduced_vector_1, reduced_vector_2)
    return r

def top_X_bottom_X_truncated_motifs_vector_pearsons_correlation(vector_1, vector_2, X):
    """TAKES BOTH THE TOP X AND THE BOTTOM NOTE NON-ZERO X VALUES OF EACH VECTOR,
    CONSTRUCTS THE UNION VECTOR FOR EACH, AND THEN RETURNS THE CORRELATION SCORE"""
    top_indices_1 = vector_1.nlargest(X).index
    top_indices_2 = vector_2.nlargest(X).index
    bottom_indices_1 = vector_1[vector_1 != 0].nsmallest(X).index
    bottom_indices_2 = vector_2[vector_2 != 0].nsmallest(X).index
    combined_indices_1 = top_indices_1.union(bottom_indices_1)
    combined_indices_2 = top_indices_2.union(bottom_indices_2)
    combined_indices = combined_indices_1.union(combined_indices_2)
    reduced_vector_1 = vector_1[combined_indices]
    reduced_vector_2 = vector_2[combined_indices]
    r = pcc(reduced_vector_1, reduced_vector_2)
    return r

def convert_dataframe_to_MEGA11_file_format(distance_df,title):
    """TAKES THE DISTANCE MATRIX (STORED AS A DATAFRAME) AND CONVERTS IT INTO A STRING IN THE EXACT .MEG FILE FORMAT WITHOUT ANY ERRORS"""
    n_rows, n_cols = distance_df.shape
    file_string = f"#mega\n!Title {title} distance matrix;\n!Format DataType=Distance;\n\n"
    animal_species_list = distance_df.columns
    for animal in animal_species_list:
        file_string += f"#{animal}\n"
    file_string += "\n\n"
    row_ind_string = "\t[0"
    for i in range(1,n_rows):
        row_ind_string+=f"\t{i}"
    file_string+= f"{row_ind_string}]\n"

    for i in range(n_rows):
        new_row = f"[{i}]\t"
        for j in range(i): 
            new_row += f"{distance_df.iloc[i,j]}\t"
        file_string+=f"{new_row}\n"
    return file_string

def write_string_to_file(string, filepath):
    """WRITES THE STRING TO THE MENTIONED FILEPATH"""
    try:
        with open(filepath, 'w') as file:
            file.write(string)
        print(f"String successfully written to {filepath}")
    except Exception as e:
        print(f"An error occurred: {e}")

def generate_variable_Nmers(N):
    """
    Generate all possible variable k-mers for a given k.
    The first and last positions are fixed (A, G, C, T), while intermediate positions can be 'N'.
    """
    if N < 3:
        raise ValueError("k must be at least 3 to allow fixed start/end positions and variable intermediates.")
    def extend_iterative_Nmer_list(kmer_list):
        """TAKES EACH KMER IN THE GIVEN LIST, ADDS AN EXTRA NUCLEOTIDE TO ITS END FOR EACH SUCH NUCLEOTIDE,
        AND THEN RETURNS THE FINAL LIST OF ALL SUCH (K+1)MERS AS A PYTHON LIST"""
        extended_Nmer_list = []
        residues = ['A','T','G','C','N']
        for kmer in kmer_list:
            for nucl in residues:
                extended_Nmer_list.append(kmer+nucl)
        return extended_Nmer_list
    nucleotides = ['A','T','G','C','N']
    output = nucleotides
    for i in range(N-3):
        output = extend_iterative_Nmer_list(output)
    kmer_list = generate_kmer_motif_list(N)
    penultimate_output = []
    for letter in ['A','T','G','C']:
        for item in output:
            penultimate_output.append(letter+item)
    ultimate_output = []
    for letter in ['A','T','G','C']:
        for item in penultimate_output:
            ultimate_output.append(item+letter)
    return [item for item in ultimate_output if item not in kmer_list]


'''
#---------------------------------------------------------------------------------------------------------

def synteny_blocks_wrt_one_gene_order(ref_gen_order, analysis_gen_order):
    synteny_blocks = []
    block = []
    i = 0
    dir = 0
    while i<len(analysis_gen_order):
        if len(block)==0:
            block.append(analysis_gen_order[i])
            i+=1
        if len(block)==1:
            j = ref_gen_order.index(block[0]) #assuming block[0] is in ref_gen_order!!!
            if ref_gen_order[j+1]==analysis_gen_order[i]:
                block.append(ref_gen_order[j+1])
                dir = 1
                i+=1
            elif ref_gen_order[j-1]==analysis_gen_order[i]:
                block.insert(0,ref_gen_order[j-1])
                dir = -1
                i+=1
            else:
                synteny_blocks.append(block)
                block = []
'''
        

def identify_synteny_blocks_between_two_gene_orders_with_same_gene_content(gene_order_1, gene_order_2):
    synteny_blocks_1 = []
    synteny_blocks_2 = []

all_species_names = list_all_species_names_from_file_path()
all_species_gene_data = {species: pd.read_csv(f"data/genbank_files/{species}/{species}_cleaned_gene_data.tsv", 
                                              delimiter="\t").squeeze() for species in all_species_names}
all_species_gene_order_lists = {species: all_species_gene_data[species]['Gene'].tolist() 
                                for species in all_species_names}

human_gene_order = all_species_gene_order_lists['Homo_sapiens']
fruitfly_gene_order = all_species_gene_order_lists['Drosophila_melanogaster']
for item in synteny_blocks_wrt_one_gene_order(human_gene_order, fruitfly_gene_order):
    print(item)
    print("")