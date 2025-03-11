# mitoevo
This repo contains all code, docs and other files that I have used so far whilst studying and analysing mitochondrial genomes for evidence and patterns behind translocations, using ideas underlying motif distributions.

1. First run an NCBI query search and download the NCBI IDs of the complete mitochondrial genomes of interest along with their scientific names. Store these in a .tsv file as two columns corresponding to ID and Name respectively.
   \\ I need to look for whether I can add a third column classfying them into phyla/similar groups
2. Then, `a1` downloads the genbank files for all of these listed IDs and saves them according to their species names in folders with the same species name.
3. We then need to clean up the gene data and rename all mitochondrial gene orders to a standard naming convention. For this purpose, `a2.0` takes in all possible `Product` to `Gene` mappings and stores it as a .tsv file.
4. Then, we need to clean up the saved mapping file *MANUALLY*, following which `a2.1` will refer to this file and create .tsv files for each genbank file such that the gene names are all filled in cleanly
   \\ The only catch here is that it won't differentiate between multiple tRNA genes in a single genome. What I could then do is to filter out my dataset in such a way that only CDS-gene orders with number of genomes > 3 will be retained, while all other anomalous gene orders will be deleted due to lack of scope of tractability.
5. Then, we can play around once we have these gene order .tsv files saved. This is what `a3.0` and `a3.1` try to do. I stil have to figure out what these files exactly produce and what are some nice ways to visualise the data that we have. 
6. Lastly, `a4` will take all of these genomes and construct whole-genome vectors for all of them, both 5-mer and 6-mer vectors. This will be then used to construct trees. Again, I'll have to choose my dataset in such a way that I can get some meaningful results and inferences out of the clusters that we see forming in the tree from the vectors versus the tree from using an MSA.

7. Then, we proceed to main the algorithm itself, which is to obtain the sequences/regions we are after. We wish to take all of the sequences in our dataset and assign distances between each of them pair-wise. For this, we need to construct our own metric, by borrowing some hints from the `Levenshtein distance` and modifying it for our own cause. However, we need to get an idea of how many CDS-orders are there in our dataset, so that we can play around with this and see how much of a 'weight' we need to assign to the CDS-gene order itself over the RNA-gene order.
8. Once we have the Levenshtein-Lokesh metric distances between all the genomes, we can then construct the trees based on the distance matrix and compare with the MSA as well as the vector alignment trees to see how the clusters have formed and whether they match or not.
9. Next, using this same distance matrix, we need to come up with some `order of evolution` for all the genomes in this dataset. *I need to figure out how to do this from the distance matrix...*
10. Between each pair of genomes, you then need to identify which chunks of genes have changed between two gene orders. *I need to come up with a clever function for this. The old stuff might not work. In fact the older function is probably incorrect and has a lot of errors in it...*
11. Then, for each pair of genomes as I iterate over the dataset, and for each translocated gene chunk, I identify the four sets of sequences using the following hyperparameters:
    1.  The window size from the 5' end (before the left side) = $100$ (prolly)
    2.  The window size from the 3' end (before the right side) = $100$ (prolly)
    I will then need to optimise over these hyperparameters later on (*I don't know what my optimisation criterion would be*) and find the best hyperparameters.
    
12. For each of these four categories of sequences obtained, I will then need to perform further analysis in the following manner. I first treat these four categories of sequences separately, then as two pairs, and finally club them together in the following analysis. For each of these sequences, I construct the 5-mer and 6-mer vectors for them and store all of this data. 
13. Once I have the vectors for each of these sequences, I can then do more analysis:
    1.  I can look for motifs that are common to all of these vectors (or at least, occur in >60% of these vectors)
    2.  I can look for specific motif sequences based on literature about existing known motifs for chromosomal translocations
    3.  I can do a `truncated-correlation` and use that maybe as my `optimisation metric` for the optimisation part
14. Further from this analysis we can compare what motifs we get and see if there are any unexpected ones that we need to focus on... but yes, until this would be a good part of work that we get done with.

:)