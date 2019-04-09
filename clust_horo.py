#!/usr/bin/env python3

#Import libraries
import argparse
import sys
import os
import Bio

#Function used to grab the flags
def get_args():
   
   #Define all the flags
   argparser = argparse.ArgumentParser()
   argparser.add_argument("-i", "--input", required = True, help = "Input file name is required! Input file must be in fasta format.")
   argparser.add_argument("-t", "--threshold", type = float, default = 0.5, help = "The minimum threshold desired when clustering the input sequences")
   argparser.add_argument("-o", "--output", type = str, help = "Output file name to which the results are written to")
   argparser.add_argument("-m", "--minseqlen", default = 0, type = int, help = "Throw out any input sequences that are below the minimum sequence length")   
   argparser.add_argument("-s", "--sort", type = int, default = 0,  help = "Select 0 to not sort the input (default), 1 to sort the input file in decreasing sequence length, 2 to sort by decreasing abundance, and 3 to sort by both")
   
   args = argparser.parse_args()

   # If an output filename is not specified
   if not args.output:

      # Assign a generic output name
      args.output = 'output_'+args.input

   return args

#Function for running cluster_fast without any sorting
def sort_default(filename, t, minseqlen):
   cmd = 'usearch -cluster_fast '+ filename + 'id ' +str(t)+ ' -clusters cluster -centroids centroids.txt -minseqlength ' +str(minseqlen)
   os.system(cmd)

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
    
    cmd = 'usearch -cluster_fast sorted_abundance.txt -id '+str(t)+' -clusters cluster -centroids centroids.txt -minseqlength ' + str(minseqlen)
    os.system(cmd)    
    #os.system('rm sorted_abundance.txt')


# Function for sorting by length
def sort_by_length(filename, t, minseqlen):
    cmd = 'usearch -cluster_fast '+filename+' -id '+str(t)+' -clusters cluster -centroids centroids.txt -sort length -minseqlength ' + str(minseqlen)
    os.system(cmd)

# Function for sorting by abundance and length
def abundance_and_length(filename, t, minseqlen):
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
    
    sort_by_length('sorted_abundance.txt', t, minseqlen)
    #os.system('rm sorted_abundance.txt')

def clust_homo():
   file = open("centroids.txt")
   data = file.read()

def main():
   args = get_args()
   
   file = open(args.input).read()
   
   filename_input = args.input
   t = args.threshold
   minseqlen = args.minseqlen

   if args.sort == 0:
      sort_default(filename_input, t, minseqlen)

   if args.sort == 1:
      sort_by_length(filename_input, t, minseqlen)
      
   elif args.sort == 2:
      sort_by_abundance(filename_input, t, minseqlen)
   
   else:
      abundance_and_length(filename_input, t, minseqlen)

   clust_homo()

if __name__ == "__main__":
   main()
