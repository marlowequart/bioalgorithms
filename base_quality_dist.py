'''
Given: FASTQ file, quality threshold q


Return: Number of positions where mean base quality falls below given threshold
'''

from Bio import SeqIO

file='/Users/Marlowe/Downloads/rosalind_bphr.txt'

records=SeqIO.parse(file,"fastq")
sequence_dict=SeqIO.to_dict(SeqIO.parse(file,"fastq"))

with open(file) as file:
	quality=int(float(file.readline()))
	
	
#create a list of a list of the phred quality scores from each record
quality_scores=[]
for record in records:
	quality_scores.append(record.letter_annotations["phred_quality"])

#print(quality_scores[0][0])
'''
find number of phred quality scores in a record
add together the nth phred quality scores in all of the records
divide that number by the number of records to give the avearge of the nth base quality score
'''
#get the number of quality scores in the records
num_quality_scores_per_record=len(quality_scores[0])
num_records=len(quality_scores)
#Find the mean for each read
	
current_read_mean2=[]
current_read_mean2.append([sum(quality_scores[y][x] for y in range (0,num_records))/
							num_records for x in range(0,num_quality_scores_per_record)])

count=0
for i in current_read_mean2[0]:
	if i < quality:
		count += 1
	
	
print(count)