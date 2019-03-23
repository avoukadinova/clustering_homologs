#!/usr/bin/env python3
import argparse
import sys

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

def main():
   args = get_args()
   
   f = open(args.input).read
   
   print("The input file name is " + args.input)
   print("The output file name is " + args.output)

if __name__ == "__main__":
   main()
