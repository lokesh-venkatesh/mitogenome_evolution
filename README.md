# Mitochondrial Genome Evolution
This repo contains all code, docs and other files that I have used so far whilst studying and analysing mitochondrial genomes for evidence and patterns behind translocations, using ideas underlying motif distributions.

A brief summary of the workflow followed across this project is as follows:
### Looking at motif distributions as a metric for reconstructing phylogeny 
We looked at whole-genome correlations as well as gene-wise motif distributions as a way to draw a 'similarity score' between organisms for means of constructing phylogeny. This was followed by a score-based comparison of the phylogenetic trees obtained, which we evaluated to look at how useful of a metric motif distributions are in reconstructing phylogeny. 

### Looking at regulatory correlations within the mitochondrial genome as well as across the mitochondrial and nuclear genomes:
We also looked at large non-geneic regions covering large swathes of an organism's mitochondrial genome, and looked at their motif distributions across various organisms, looking for trends and correlations. This was followed by looking at correlations with these same motif distribution vectors with regularly spaced intervals of genomic regions in the nuclear chromosomes for a select few organisms, including Homo sapiens.

### Identifying sequences that correspond to 'chopping' due to translocations over the course of mitochondrial evolution:
The following was the workflow we used to be able to isolate a number of sequences corresponding to translocation 'chop' sites:
1. Query every mitochondrial genome available on the NCBI database
2. Clean up the data and obtain 'sequences' of each genome's gene orders
3. Perform an MSA using the gene orders as the characters
4. Within each cluster of genomes after performing the MSA, perform a sequence-based MSA
5. This is then followed by a whole order-rearrangement of the 15k genomes queried
6. Then, we identify pairs of genomes where the gene orders have changed. This indicates those pairs of genomes where a translocation has taken place.
7. From each of these tuples of genomes, we identify the synteny blocks that have undergone a rearrangement, and then isolate the four types of sequences related to the rearrangement of one synteny block
8. From each of these four sequences obtained, we then compile them for all such tuples obtained, and finally analyse these sequences for common motifs and trends in their motif distributions for any 'signal' that a hypothetical translocation-causing 'machinery' would identify and execute.
