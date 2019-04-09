"""
Hierarchical Agglomerative Clustering (Centroid Linkage) Algorithm for USEARCH

Time complexity: O(n^2log(n))
Hierarchical finds non-overlapping clusters
One of the methods used in this script is agglomerative
This merges nodes with a bottom-up approach
By finding the least amount of distance between centroids.

There are four different methods in hierarchical agglomerative clustering
One of them is centroid linkage.
Centroid method or linkage is the distance between 2 mean vectors of the distance of 2 clusters
In this script, two clusters are combined to output the smallest centroid distance.

measurement of d12 between cluster 1 and cluster 2:
d12 = d(x,y)


TO DO/ADD
	+add .csv output
 		import csv, sys 
		
			with open('clusters.csv', 'w') as f:
			out = csv.writer(f)
			out.outrow(zip(row))
""""


import itertools
from numpy import mean
from Bio import SeqIO
from Bio.Align import AlignInfo
import argparse
        
# input sequences
class Clust(object):

	clust = [] #clusters
	s = [] #sequences
	
	def __init__(self, s):
		self.s = s

# hierarchical clustering class, inherits from Clust class
class HierClust(Clust):

	#define distance matrix
	dist = []
	
# Use hamming distance to leave deletion types
# This counts the number of differences between the homologs (equal length strings: s1, s2)
	def hamm(s1, s2):
		diff = 0 #difference
        len = 0 #length
        for ch1, ch2 in zip(s1, s2):
                if ch1 != ch2:
                        diffs += 1
                        len += 1
        return diff
        
# alternative using imap - fast/less memory
	# def hamming_distance(s1, s2):
		# assert len(s1) == len(s2)
		# ne = operator.ne
		# return sum(imap(ne, s2, s2))

	# iterates clusters across each cluster
	def __init__(self, cent, thresh):
		Clust.__init__(self, seqs)
		self.cent = cent
		self.thresh = thresh
		self.seqs = list()
		if self.cent.name != 'consensus':
			self.seqs.append(cent)
			
	def __repr__(self):
		return "\nCluster: %s Centroid: %s\n" % (self.size(), self.cent.seq)
		
	# Hamming distance that connects the center and sequences
	def dist(self, doc):
		return hamm(self.cent.seq, doc.seq)
		
	# Add a cluster member to this specific cluster	
	def add(self, doc):
		if doc.name != 'consensus':
			self.seqs.append(doc)
			
	# Two clusters will merge
	def two_merge(self, clust):
		self.seqs += clust.seqs
		
	# This averages the distance to the other cluster members	
	def dist_a(self, doc):
		return mean(hamm(r.seq, doc.seq) for r in self.seqs)
		
	# This allocates a new centroid to a cluster
	def center(self):
		self.cent = self.consensus()
		return self.cent
		
	# this counts the number of cluster members	
	def size(self):
		return len(self.seqs)
	
	def consensus(self): 
	# Return consensus sequence
		a = Align.MultipleSeqAlignment(self.seqs)
		i = AlignInfo.SummaryInfo(a)
		s = info.dumb_consensus(self, threshold=.7, ambiguous="X", consensus_alpha=None, require_multiple=0)
			# s outputs a fast consensus sequence of the alignment
			# this will go back tot he sequence residue
			# and calculate the number of each type of residue in all of the sequence alignment
		
		return SeqRecord.SeqRecord(seq, name='consensus', id='consensus')

	# iterat variable will iterate on all clusters
	# this yields each other the cluster centroids
	# then yields the next most central part from every cluster from its cluster sequence
	# until it has gone through all of the cluster sequences
	def iterat(self):
		self.iter()
    		yield self.cent
    		for_, doc in sorted((self.dist(r), r) for r in self.seqs):
    			yield doc
    
    # this will take the input		
    def handler(args):
    	seqrecords = SeqIO.parse(args.alignment, 'fasta')
    	seqrecords = (x for _, x in sorted((-len(sr.seq.ungap('-')), sr) for sr in seqrecords))
    	
    	# cluster
    	clust = Clust(seqrecords, args.threshold, args.thresh)
    	print "Clustering complete!"
    	if args.min_per_clust:
    		clust.merge_small_clusters(args.min_per_clust)
    		print "Merge small clusters: DONE"
    		
    	#output	-csv file add?
    	cluster.write(args.output)
    	


if __name__ == '__main__':
    main(get_args())
