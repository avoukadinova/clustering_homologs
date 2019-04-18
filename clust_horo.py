#!/usr/bin/env python3

# Import libraries
import argparse
import sys
import subprocess
import os
from Bio import SeqIO
from itertools import combinations
from Bio import pairwise2
import time
from itertools import combinations

# Used as the intitial time for clust_horo run time
t0 = time.time()

# Function used to grab the flags
def get_args():
   
    # Define all the flags
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-i", "--input", required = True, help = "Input file name is required! Input file must be in fasta format.")
    argparser.add_argument("-t", "--threshold", type = float, default = 0.5, help = "The minimum threshold desired when clustering the input sequences")
    argparser.add_argument("-o", "--output", type = str, help = "Output file name to which the results are written to. There is no need to provide an extension; this is simply the name you want on your clusters. An example of a common output name is 'cluster_genomename'")
    argparser.add_argument("-m", "--minseqlen", default = 0, type = int, help = "Throw out any input sequences that are below the minimum sequence length")   
    argparser.add_argument("-s", "--sort", type = int, default = 0,  help = "Select 0 to not sort the input (default), 1 to sort the input file in decreasing sequence length, 2 to sort by decreasing abundance, and 3 to sort by both")
    argparser.add_argument("-p", "--processors", type = int, default = 1, help = "Input the total number of processors you want to use.")

   
    args = argparser.parse_args()

    # If an output filename is not specified
    if not args.output:
        # Assign a generic output name
        args.output = 'output_'+args.input

    return args

# Function for running cluster_fast without any sorting
def sort_default(filename, t, minseqlen, p):
    
    # Used as the initial time for USEARCH run time
    t2 = time.time()

    # Run USEARCH
    cmd = 'usearch -cluster_fast '+ filename + ' -id ' +str(t)+ ' -clusters cluster -centroids centroids.txt -minseqlength ' +str(minseqlen) + ' -threads ' + str(p)
    os.system(cmd)

    # Print the total USEARCH run time
    t3 = time.time()
    print("USEARCH run time: " + str(t3-t2) + " seconds")

# Function for sorting by abundance
def sort_by_abundance(filename, t, minseqlen, p):

    cmd = 'usearch -cluster_fast '+filename+' -id 1.0 -clusters cluster -centroids centroids.txt -threads ' + str(p)
    os.system(cmd)
    
    handle=open('centroids.txt').read()
    centroid_count = handle.count('>')
    
    clusters=[]
    for x in range(centroid_count):
        cluster = open('cluster'+str(x)).read()
        freq = cluster.count('>')
        clusters.append([freq, 'cluster'+str(x)])
    
    clusters.sort(key = lambda tup: tup[0], reverse=True)
    
    handle = open('sorted_abundance.txt','a')
    
    for cluster in clusters:
        name = cluster[1]
        temp = open(name).read()
        handle.write(temp)
    
    os.system('rm cluster*')
    os.system('rm centroids.txt')
    handle.close()

    t2 = time.time()
    
    cmd = 'usearch -cluster_fast sorted_abundance.txt -id '+str(t)+' -clusters cluster -centroids centroids.txt -minseqlength ' + str(minseqlen) + ' -threads ' + str(p)
    os.system(cmd)    

    t3 = time.time()

    print("USEARCH run time: " + str(t3-t2) + " seconds.")

# Function for sorting by length
def sort_by_length(filename, t, minseqlen, p):
   
    t2 = time.time()

    cmd = 'usearch -cluster_fast '+filename+' -id '+str(t)+' -clusters cluster -centroids centroids.txt -sort length -minseqlength ' + str(minseqlen) + ' -threads ' + str(p)
    os.system(cmd)
    
    t3 = time.time()
    
    print("USEARCH run time: " + str(t3-t2) + " seconds.")

# Function for sorting by abundance and length
def abundance_and_length(filename, t, minseqlen, p):

    cmd = 'usearch -cluster_fast '+filename+' -id 1.0 -clusters cluster -centroids centroids.txt -minseqlength ' + str(minseqlen) + ' -threads ' + str(p)
    os.system(cmd)
    
    handle=open('centroids.txt').read()
    centroid_count = handle.count('>')
    
    clusters=[]
    for x in range(centroid_count):
        cluster = open('cluster'+str(x)).read()
        freq = cluster.count('>')
        clusters.append([freq, 'cluster'+str(x)])
    
    clusters.sort(key = lambda tup: tup[0], reverse=True)
    
    handle = open('sorted_abundance.txt','a')
    
    for cluster in clusters:
        name = cluster[1]
        temp = open(name).read()        
        handle.write(temp)
    
    os.system('rm cluster*')
    os.system('rm centroids.txt')
    handle.close()
    
    sort_by_length('sorted_abundance.txt', t, minseqlen, p)

# Function for moving the USEARCH output to a new USEARCH directory
def move_files():
    
    # Check if there is a results directory and if not, make one
    if os.path.isdir("./USEARCH_results") != 1:
        os.system("mkdir ./USEARCH_results")
   
    # Move the USEARCH output to the new USEARCH directory
    os.system("mv ./cluster1* ./USEARCH_results/")
    os.system("mv ./cluster2* ./USEARCH_results/")
    os.system("mv ./cluster3* ./USEARCH_results/")
    os.system("mv ./cluster4* ./USEARCH_results/")
    os.system("mv ./cluster5* ./USEARCH_results/")
    os.system("mv ./cluster* ./USEARCH_results/")
    #os.system("mv ./centroids.txt ./USEARCH_results/")

    # Remove unwanted files
    if os.path.exists("./sorted_abundance.txt") == 1:
        os.system("rm sorted_abundance.txt")

                 
def cluster_compare(i, j):
    cluster1 = []
    cluster2 = []

    for record in SeqIO.parse("./USEARCH_results/cluster" + str(i), "fasta"):
        cluster1.append(record.seq)

    for record in SeqIO.parse("./USEARCH_results/cluster" + str(j), "fasta"):
        cluster2.append(record.seq)

    for sequence1 in cluster1:
        for sequence2 in cluster2:
            sim_score = similarity(sequence1, sequence2)          
            if sim_score >= t:
                continue
            else:
                return False

    return True

# similarity measure using pairwise global alignment (BLAST definition)
def similarity(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    alignments = pairwise2.align.globalxx(seq1, seq2)
    score = alignments[0][2]
    columns = len(alignments[0][0])

    return score / columns

def clust_horo(t):

    # if os.path.isdir("./clust_horo/") != 1:
    #     os.system("mkdir ./clust_horo/")   

    centroids = []
    # num_centroids = int(subprocess.getoutput('ls | grep "cluster" | wc -l')) gotta figure out how to get correct pwd for this
    num_centroids = open('centroids.txt').read().count('>')
    #print(num_centroids)

    for i in range(num_centroids):
        centroid = next(SeqIO.parse('./USEARCH_results/cluster'+str(i), 'fasta'))
        centroids.append(centroid.seq)             

    cluster_num = 0
 
    with open('centroids_sorted.txt', 'w') as f:
        for item in centroids:
             f.write(">cluster"+ str(cluster_num) + "\n") 
             f.write(str(item) + "\n")
             cluster_num += 1

    f.close()
 
    cmd = 'usearch -cluster_fast centroids_sorted.txt -id '+str(t)+' -centroids centroids_clustered.txt'
    os.system(cmd)

    num_new_centroids = open('centroids_clustered.txt').read().count('>')
    #print(num_new_centroids)    
    print("clust_horo found " + str(num_centroids - num_new_centroids) + " homologous clusters!")

    for i in range(num_new_centroids):
         if (">cluster" + str(i) + '\n') not in open("centroids_clustered.txt").read():
              print("cluster" + str(i))
    #for centroid1, centroid2 in combinations(centroids, 2):

    #    sim_score = similarity(str(centroid1), str(centroid2))

    #    if sim_score < t:
    #        continue
    #    else:
    #        print("They might be the same\n")

    #        i = centroids.index(centroid1)
    #        j = centroids.index(centroid2)
            # do this here: 'I’d keep a list of all the skipped #s to avoid going through the clusters again at the end' 
    #        same = cluster_compare(i, j)

    #        if same == False:
    #            print("cluster" + str(i) + ' and cluster' + str(j) + ' should be one cluster')
                # cmd = 'cat cluster' + str(i) + ' cluster' + str(j) + ' > cluster' + str(i)
                # os.system(cmd)

                # cmd = 'rm cluster' + str(j)
                # os.system(cmd)

                # cmd = 'mv cluster' + str(len(centroids)) + ' cluster' + str(j)
                # os.system(cmd)

                # what about finding new centroid and removing old ones?

    #       else:
    #           continue

def main():
    
    # Get the flags
    args = get_args()
   
    # Grab the input filename, threshold, and minseqlen
    filename_input = args.input
    t = args.threshold
    minseqlen = args.minseqlen
    p = args.processors
                 
    # Run the desired sorting algorithm and USEARCH

    if args.sort == 0:
        sort_default(filename_input, t, minseqlen, p)

    elif args.sort == 1:
        sort_by_length(filename_input, t, minseqlen, p)
      
    elif args.sort == 2:
        sort_by_abundance(filename_input, t, minseqlen, p)
   
    else:
        abundance_and_length(filename_input, t, minseqlen, p)

    move_files()
    clust_horo(t)

    # Used as the final time for the clust_horo total run time
    t1 = time.time()
    print("clust_horo total run time: " + str(t1-t0) + " seconds.\n")

if __name__ == "__main__":
    main()
