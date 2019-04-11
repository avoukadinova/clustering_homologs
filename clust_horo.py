#!/usr/bin/env python3

#Import libraries
import argparse
import sys
import os
from Bio import SeqIO
import time

#Used as the intitial time for clust_horo run time
t0 = time.time()

#Function used to grab the flags
def get_args():
   
   #Define all the flags
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

#Function for running cluster_fast without any sorting
def sort_default(filename, t, minseqlen):

   #Used as the initial time for USEARCH run time
   t2 = time.time()
   
   #Running USEARCH
   cmd = 'usearch -cluster_fast '+ filename + ' -id ' +str(t)+ ' -clusters cluster -centroids centroids.txt -minseqlength ' +str(minseqlen)
   os.system(cmd)

   #Print the total USEARCH run time
   t3 = time.time()
   print("USEARCH run time: " + str(t3-t2) + " seconds")

#Function for sorting by abundance
def sort_by_abundance(filename, t, minseqlen):
    cmd = 'usearch -cluster_fast '+filename+' -id 1.0 -clusters cluster -centroids centroids.txt'
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
    
    cmd = 'usearch -cluster_fast sorted_abundance.txt -id '+str(t)+' -clusters cluster -centroids centroids.txt -minseqlength ' + str(minseqlen)
    os.system(cmd)    

    t3 = time.time()

    print("USEARCH run time: " + str(t3-t2) + " seconds.")

# Function for sorting by length
def sort_by_length(filename, t, minseqlen):
    t2 = time.time()

    cmd = 'usearch -cluster_fast '+filename+' -id '+str(t)+' -clusters cluster -centroids centroids.txt -sort length -minseqlength ' + str(minseqlen)
    os.system(cmd)
    
    t3 = time.time()
    
    print("USEARCH run time: " + str(t3-t2) + " seconds.")

# Function for sorting by abundance and length
def abundance_and_length(filename, t, minseqlen):
    cmd = 'usearch -cluster_fast '+filename+' -id 1.0 -clusters cluster -centroids centroids.txt -minseqlength ' + str(minseqlen)
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
    
    sort_by_length('sorted_abundance.txt', t, minseqlen)

#Function for moving the USEARCH output to a new USEARCH directory
def move_files():

   #Check if there is a USEARCH directory and if not, make one
   if os.path.isdir("./USEARCH") != 1:
       os.system("mkdir ./USEARCH")
   
   #Move the USEARCH output to the new USEARCH directory
   os.system("mv ./cluster* ./USEARCH/")
   os.system("mv ./centroids.txt ./USEARCH/")

   #Remove unwanted files
   if os.path.exists("./sorted_abundance.txt") == 1:
       os.system("rm sorted_abundance.txt")

#Function for splitting a sequence into length 'size'
def split_seq(seq,size):
  
    return [seq[i:i+size] for i in range(0, len(seq), size)]

#Function for finding the Jaccardian distance
def jaccard(s1, s2):
    
   #print("TO DO! Function for findind Jaccardian distance! Currently does not work.") 
   return 0

#Function for reclustering the USEARCH clusters
def clust_horo():

   if os.path.isdir("./clust_horo/") != 1:
       os.system("mkdir ./clust_horo/")   
   
   centroids = []
   centroids_copy = []

   for record in SeqIO.parse("./USEARCH/centroids.txt", "fasta"):

       centroids.append(record)
       centroids_copy.append(record)              

   for centroid1 in centroids:
       for centroid2 in centroids_copy:

           s1 = split_seq(centroid1, 3)
           s2 = split_seq(centroid2, 3)

           sim_score = jaccard(s1, s2)          
           #print(sim_score)

def main():

   #Get the flags
   args = get_args()
   
   #Grab the input filename, threshold, and minseqlen
   filename_input = args.input
   t = args.threshold
   minseqlen = args.minseqlen

   #Run the desired sorting algorithm and USEARCH
   if args.sort == 0:
      sort_default(filename_input, t, minseqlen)

   elif args.sort == 1:
      sort_by_length(filename_input, t, minseqlen)
      
   elif args.sort == 2:
      sort_by_abundance(filename_input, t, minseqlen)
   
   else:
      abundance_and_length(filename_input, t, minseqlen)

   move_files()
   clust_horo()

   #Used as the final time for the clust_horo total run time
   t1 = time.time()
   print("clust_horo total run time: " + str(t1-t0) + " seconds.\n")

if __name__ == "__main__":
   main()
