#!/usr/bin/env python3

# Import libraries
import argparse
import sys
import subprocess
import os
from Bio import SeqIO
import time

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

def clust_horo(t):   

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
 
    cmd = 'usearch -cluster_fast centroids_sorted.txt -id '+str(t)+' -centroids centroids_clustered.txt -clusters new'
    os.system(cmd)

    num_new_centroids = open('centroids_clustered.txt').read().count('>')
    #print(num_new_centroids)    
    print("clust_horo found " + str(num_centroids - num_new_centroids) + " homologous clusters!\n")
   
    f_out = open("homo_out.txt", "w")
    for i in range(num_new_centroids):

         filename = "new" + str(i)
         num_seq = open(filename).read().count('>')

         if(num_seq > 1):
            cluster = open(filename).read()
            sequences = SeqIO.parse(filename, 'fasta')

            clusters = []
            for seq in sequences:
               clust_name = seq.id
               clust_num = clust_name.split("cluster")[1]
               clusters.append(clust_num)
               f_out.write(seq.id)
               f_out.write(" ")

            f_out.write("\n")
   
    f_out.close()

    os.system("rm new1*")
    os.system("rm new2*")
    os.system("rm new3*")
    os.system("rm new4*")
    os.system("rm new5*")
    os.system("rm new6*")
    os.system("rm new7*")
    os.system("rm new*")
    
    os.system("rm centroid*")

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
