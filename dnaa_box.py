'''
DNAA Boxes

The aim of this script is to find the origin of replication for a given genome

Given:
string = a region of a genome where the supposed origin of replication is located as
suggested by finding the minimum of the skew for the genome.
k = an estimated size of dnaA string for genome being tested.
d = an estimated number of mismatches to find the neighborhoods of the dnaA boxes

Return:
a list of kmers of length k that are the potential dnaA boxes in the given region

'''

import time
from collections import defaultdict
import itertools


def get_data():
	#Load genome Data
	with open('...test_genome.txt') as file:
		title=file.readline().replace('\n','')
		seq_list=file.readlines()
		file.close()
	
	string=''.join(seq_list).replace('\n','')
	return string 


def rev_comp(text):
	#output the reverse complement of the given string
	output_forward=[]
	output=[]
	output_list=[]
	for i in text:
		if i=='C':
			output_forward.append('G')
		if i=='T':
			output_forward.append('A')
		if i=='G':
			output_forward.append('C')
		if i=='A':
			output_forward.append('T')

	for j in range(1,len(text)+1):
		output_list.append(output_forward[-j])
		
	output=''.join(let for let in output_list)
	
	return output


def HammDist(seq1,seq2):
	#The hamming distance between two strings is the
	#number of mismatches between the two strings
	count=0
	for i in range(len(seq1)):
		if seq1[i]!=seq2[i]:
			count=count+1
			
	return count

def Neighbors(pattern,d):
	'''
	finds the neighbors with hamm dist < d recursively
	'''
	if d==0:
		return [pattern]
	if len(pattern)==1:
		return ['A','C','G','T']
	neighborhood=[]
	SuffixNeighbors=Neighbors(pattern[1:],d)
	for kmer in SuffixNeighbors:
		if HammDist(kmer,pattern[1:]) < d:
			for x in ['A','C','G','T']:
				neighborhood.append(x+kmer)
		else:
			neighborhood.append(pattern[0]+kmer)
	return neighborhood


def find_kmer(string,k,d):
	'''
	In an effort to speed up the process of finding the most frequent kmers,
	With this version, we will generate the neighborhood for each kmer in sequence,
	reducing the number of computations significantly.
	
	This function outputs a list of the kmers with the maximum number of neighbors
	that have a hamming distance < d
	'''
	#First create a list of all the kmers of length k in string and rev complement of string
	kmer_dict=[]
	for i in range(len(string)-k+1):
		kmer_dict.append(string[i:i+k])
		

		
	#Create an array of all possible kmers of length k, store in dictionary
	arrays=list(itertools.product(['A','C','G','T'],repeat=k))
	outs={}
	outs=defaultdict(lambda:0,outs)	
	for item in arrays:
		string=''.join(let for let in item)
		outs[string]


	#Next find the neighborhood with hamm dist<d for each kmer in string
	#save the kmers that meet hamm dist to the dictionary	
	for kmer in kmer_dict:
		#add rev comp of kmer here to get neighbors of both.
		nbrs=Neighbors(kmer,d)+Neighbors(rev_comp(kmer),d)
		for item in nbrs:
			outs[item] += 1
			
	#find maximum value of kmers with hamm dist<d
	maxi=0
	for key in outs.keys():
		if outs[key]>maxi:
			maxi=outs[key]

	#create list of kmers with that maximum value
	kmers_out=[]
	for key in outs.keys():
		if outs[key]==maxi:
			kmers_out.append(key)
			
	return kmers_out
		
	
	
def save_output_str(output_list):
	''' Convert a list of strings to a string with each item in list seperated by newline
		save string to text file
	'''
	string=' '.join(str for str in output_list)
	with open('dna_boxes_output.txt','w') as newfile:
		newfile.write(string)
		newfile.close()	
	
def main():
	start_time = time.time()
	string=get_data()
	
	#test a subset of the genome string with a given window size
	#minimum skew for E. coli is found at position 3923620
	start=3923620
	window=1000
	test_str=string[start-window//2:start+window//2]
	
	#find_kmer returns a list of the kmers of length k with the
	#maximum number of neighbors and that have a hamming distance < d
	#bacterial DnaA boxes are usually 9 nucleotides long, so use k=9 for bacteria
	kmers=find_kmer(test_str,9,1)
	print(kmers)
	save_output_str(kmers)
	print('%f seconds' % (time.time() - start_time))

main()
