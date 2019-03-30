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

def randomize_input(input_data):
   print("TODO! Function for randomizing the input file!")
   print("TODO! Import Biopython and parse fasta input")

def run_usearch(t, filename):
   t = str(t)
   os.system(str("usearch -cluster_fast " + filename + " -centroids centroids.txt -clusters cluster -id " + t))

def clust_homo():
   print("TODO! Re-clustering using our own clustering algorithm.")

def main():
   args = get_args()
   
   file = open(args.input).read()

   filename_input = randomize_input(file)
   filename_input = "metap2_homologs.txt"

   run_usearch(args.threshold, filename_input)

   clust_homo()

if __name__ == "__main__":
   main()
