
# Clustering Homologs
## By Adelina Voukadinova, Brian Ho, Athina Gerodias
### Overview
  When given a plethora of gene sequences (to be exact, over 100,000) from various organisms, the mere idea of determining which sequences are most similar can make one’s head (and computer) crash. Comparative biology and the ability to recognize homologous (similar in position, structure, and evolutionary origin but not necessarily in function) genes (i.e., homologs) is crucial to genome annotation, protein structure, and phylogenetics. Due to the high demand of a fast and accurate algorithm for clustering homologs, our team is going to tackle the issue by using a generally fast algorithm to form initial homolog clusters, and then further rearrange those clusters to improve accuracy.
 
 Overall, there are three classes of clustering algorithms: hierarchical, greedy heuristic, and Bayesian. Here, we only focus on hierarchical and heuristic clustering. Hierarchical clustering such as ESPRIT-Tree and BLAST is computationally more demanding than heuristic models, but is generally more accurate; it computes pairwise distances by doing a local BLAST on all the sequences and then grouping them together by a threshold of similarity. On the other hand, heuristic models such as USEARCH and UCLUST are much faster than ESPRIT-Tree, but at the cost of more space and memory. The USEARCH algorithm computes pairwise similarity and clusters sequences together by using a greedy global-alignment algorithm. 
 
 Disadvantages of USEARCH
    As we know, the main disadvantage to USEARCH is that it can sometimes inaccurately cluster sequences. This happens for two reasons:
    
1. The way USEARCH clusters homologs is by selecting the first sequence as the initial sequence to compare the rest of the sequences to. This is a common approach when clustering, but it is biased because it is considering the first sequence to be the BEST sequence to compare the rest of the sequences to when there might be a better sequence in the dataset. 
    
2. The user can specify a minimum threshold that all the sequences must meet in order to be clustered together. USEARCH will always choose to make two clusters that meet the minimum threshold even if the sequences in both clusters may be homologs. Therefore, this splits up homologs into two clusters when they should be in one.

An additional problem we noted when speaking with our project leader is that let’s say we already used USEARCH to cluster our sequences. Now, if we want to add a new sequence to the cluster, it will always choose one of the already existing clusters to add the new sequence to even if it is not considered a homolog to the rest of the sequences. Our understanding of this is if we have clusters with consensus sequences of length 30 and we try to add a sequence of length 3000, it will not make a new cluster for the sequence; it will instead add it to the best existing cluster even if it is not a good fit.


### Proposed solution
- USEARCH Parameters

  Our Python wrapper will feature two major components: the UCLUST initial clustering and the re-clustering algorithm that we will write. The first step is we need to figure out the fastest clustering method using USEARCH. Because we are only using USEARCH for its speed, we can afford to have a less accurate clustering than a command that is more accurate. Additionally, we need to finalize the parameters of our USEARCH command. The most important of these issues is the threshold value (which we are thinking of having be around 80%).
    
- Python Wrapper

  Once the USEARCH performs the initial clustering, it is our team’s duty to write an algorithm to recluster the output. The first step is identifying which sequences are placed in the wrong cluster and separating them from that cluster. Once all the misidentified sequences are separated from their clusters, we will use a clustering algorithm to either place those clusters in the right cluster or make a new cluster. Once we do this, we will allow for an option where the user can input a new sequence to an existing cluster and have it cluster it properly.

### Design
In general, clustering algorithms fall into two categories: 1. fast but inaccurate 2. slow but accurate. Our goal is to design a clustering algorithm that is both fast and accurate. First, our algorithm must handle large datasets of 80,000 to 200,000 gene sequences in fasta format. To retain the “fast” component of USEARCH, our solution will first cluster the sequences using USEARCH or UCLUST. From there, the output will include at least two clusters (but probably much more if we select an effective threshold) along with a consensus sequence for each cluster. This is where we come in. Our team will design a Python wrapper to take the USEARCH output and rearrange the sequences (and possibly form new clusters) to produce a more accurate homolog clustering. Essentially, we are using USEARCH’s output as a seed (or initial clustering) for our algorithm input. This should allow for more accurate clustering while maintaining speed. 

First, we had to decide if we wanted to rewrite a clustering algorithm or use an already existing software to create our initial cluster. Due to time constraints, we decided that it would be easier to use UCLUST to define our initial cluster and then cluster for a more accurate output. Additionally, we found it would be too difficult to reverse engineer UCLUST because Robert Edgar won’t tell the world how it works.
    Second, once we settled on using UCLUST to initially cluster our sequences, we have to decide what threshold to use for the initial cluster. This is crucial to getting the “ideal” number of clusters so that there aren’t too many or too little clusters. Reference Figure 3 for more information.
Finally, we want to make sure our algorithm can properly add new sequences to an already existing cluster since that appears to be one of the main issues with UCLUST. We are hoping to determine a fast and accurate way to do this. Additionally, our algorithm will be able to rearrange UCLUST output to make it more accurate and get rid of errors.


### Milestones
Week 1
- [x] Research USEARCH and BLAST and how they work (read the two papers given to us)
- [x] Research different clustering algorithms and decide which is best for our project
- [x] Make sure the team is comfortable with GitHub
- [x] Play around with USEARCH and make sample data

Week 2
- [x] Finalize the project and select clustering algorithms
- [x] Come up with an algorithm name
- [x] Start programming the skeleton of the project 
- [x] Finish presentation #1
- [ ] Finish GitHub repo #1

Week 3
- [ ] Have a working prototype for the project (does not have to look good)
- [ ] Start writing the paper

Week 4
- [ ] Finish GitHub repo #2
- [ ] Finish the Abstract and Introduction (paper)

Week 5
- [ ] Finish a rough draft of the project
- [ ] Finish GitHub repo #3
- [ ] Finish the Results and Conclusion (paper)

Week 6
- [ ] Finalize the project (clean code, comments, README file, algorithm handles edge cases)
- [ ] Finalize the paper

Week 7
- [ ] Finish final presentation
