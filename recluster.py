#!/usr/bin/env python3
#
from Bio import SeqIO
#from itertools import combinations
from Bio import pairwise2
import time
import csv
import os

# similarity measure using pairwise global alignment (BLAST definition)
def similarity(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    alignments = pairwise2.align.globalxx(seq1, seq2)
    score = alignments[0][2]
    columns = len(alignments[0][0])

    return score / columns

def cluster_compare(sequences, centroids):
	for centroid in centroids:
		for sequence in sequences:

			sim_score = similarity(centroid, sequence)

			if sim_score < t:
				return False
			else:
				continue

		return centroids.index(centroid)


with open('homo_out.txt') as handle:
    reader = csv.reader(handle, delimiter = ' ')
    output = list(reader)


for line in output:

	sequences = []
	centroids = []

	for cluster in line:

		records = list(SeqIO.parse('./USEARCH_results/' + cluster, 'fasta'))
		for record in records:
			sequences.append(record.seq)
		
		centroids.append(records[0].seq)

	
	result = cluster_compare(sequences, centroids)
	if result != False:
		print('These clusters should be combined.')
		print('The centroid from ' + centroids[result] + ' is the new centroid.')
		cluster1 = './USEARCH_results/' + line[0]
		cluster2 = './USEARCH_results/' + line[1]
		cmd = 'grep ' + cluster2 + '>> ' + cluster1
		os.system(cmd)
		os.system('rm ' + cluster2)

		centroids_file = SeqIO.parse('centroids_sorted.txt', 'fasta')
		file = open('centroids_sorted.txt', 'w')

		if result == 0:
			result = 1
		else:
			result = 0

		for centroid in centroids_file:
			if centroid.id = centroids[result]:
				continue:
			else:
				file.write('>' + centroid.id + '\n')
				file.write(str(centroid.seq) + '\n')

		file.close()

	else:
		#print(' and '.join(line) + 'will remain separate.')
		continue







			


