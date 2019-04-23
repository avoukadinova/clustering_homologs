#!/usr/bin/env python3
#
from Bio import SeqIO
from Bio import pairwise2
import time
import csv
import os

with open('homo_out.txt') as handle:
   reader = csv.reader(handle, delimiter = ' ')
   output = list(reader)

for line in output:

   cluster1 = './USEARCH_results/' + line[0]
   cluster2 = './USEARCH_results/' + line[1]
   cmd = 'cat ' + cluster2 + '>> ' + cluster1
   os.system(cmd)
   os.system('rm ' + cluster2)

