#!/usr/bin/env python3
import argparse
import sys
import os
import Bio

def get_args():
   argparser = argparse.ArgumentParser()
   argparser.add_argument("-i", "--input", required = True, help = "Input file name is required!")
   argparser.add_argument("-t", "--threshold", type = float, default = 0.5, help = "The minimum threshold desired for clustering the input sequences")
   argparser.add_argument("-o", "--output", type = str, help = "Output file name to which the results are written to")
   argparser.add_argument("-msl", "--minseqlen", type = int, help = "Throw out any input sequences that are below the minimum sequence length")   
 
   args = argparser.parse_args()

   if not args.output:
      args.output = 'output_'+args.input

   return args

# to sort by abundance
def sort_by_abundance(filename, t):
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
    
    cmd = 'usearch -cluster_fast sorted_abundance.txt -id '+str(t)+' -clusters cluster -centroids centroids.txt'
    os.system(cmd)    

#to sort by length
def sort_by_length(filename, t):
    cmd = 'usearch -cluster_fast '+filename+' -id '+str(t)+' -clusters cluster -centroids centroids.txt -sort length'
    os.system(cmd)

#to sort by abundance and length
def abundance_and_length(filename, t):
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
    
    sort_by_length('sorted_abundance.txt', t)

def clust_homo():
   print("TODO! Re-clustering using our own clustering algorithm.")

def main():
   args = get_args()
   
   file = open(args.input).read()
   


   clust_homo()

if __name__ == "__main__":
   main()
