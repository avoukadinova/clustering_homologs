#hierarchical/centroid clustering algorithm for USEARCH

import itertools
from numpy import mean
from Bio import SeqIO
from Bio.Align import AlignInfo
import argparse
        
# iterates clusters across each cluster
class Clust(object):

	clf = [] #cluster files
	clust = [] #clusters
	s = [] #sequences
	
	def __init__(self, s):
		self.s = s

class HierClust(Clust):

	#define distance matrix
	dist = []
	
# Use hamming distance to leave deletion types
# This counts the number of differences between the homologs (equal length strings s1, s2)
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

	def __init__(self, cent, thresh):
		Clust.__init__(self, seqs)
		self.cent = cent
		self.thresh = thresh
		self.seqs = list()
		if self.cent.name != 'consensus':
			self.seqs.append(cent)
			
	def __repr__(self):
		return "\nCluster: %s Centroid: %s\n" % (self.size(), self.cent.seq)
		
	def dist(self, doc):
		return hamm(self.cent.seq, doc.seq)
		
	def add(self, doc):
		if doc.name != 'consensus':
			self.seqs.append(doc)
			
	def two_merge(self, clust):
		self.seqs += clust.seqs
		
	def dist_a(self, doc):
		return mean(hamm(r.seq, doc.seq) for r in self.seqs)
		
	def center(self):
		self.cent = self.consensus()
		return self.cent
		
	def size(self):
		return len(self.seqs)
	
	def consensus(self): 
	#Return consensus sequence
		a = Align.MultipleSeqAlignment(self.seqs)
		i = AlignInfo.SummaryInfo(a)
		s = info.dumb_consensus(self, threshold=.7, ambiguous="X", consensus_alpha=None, require_multiple=0)
        return SeqRecord.SeqRecord(seq, name='consensus', id='consensus')
        
def iterat(self):
	self.iter()
    	yield self.cent
    	for_, doc in sorted((self.dist(r), r) for r in self.seqs):
    		yield doc
    		
    def handler(args):
    	seqrecords = SeqIO.parse(args.alignment, 'fasta')
    	seqrecords = (x for _, x in sorted((-len(sr.seq.ungap('-')), sr) for sr in seqrecords))
    	
    	clust = Clust(seqrecords, args.threshold, args.thresh)
    	print "Clustering complete!"
    	if args.min_per_clust:
    		clust.merge_small_clusters(args.min_per_clust)
    		print "Merge small clusters: DONE"
    		
    	#output	
    	cluster.write(args.output)


if __name__ == '__main__':
    main(get_args())