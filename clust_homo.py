#!/usr/bin/env python3
import argparse

def get_args():
   argparser = argparse.ArgumentParser()
   argparser.add_argument("-i", "--input", required = True, help = "Input file name is required!")
   argparser.add_argument("-t", "--threshold", type = float, default = 0.5, help = "Minimum threshold")
   return argparser.parse_args()

def main():
   args = get_args()
   
   print("The input file name is " + args.input)

if __name__ == "__main__":
   main()
