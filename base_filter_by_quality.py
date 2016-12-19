'''
Given: FASTQ file, quality cut-off value q, Phred33 quality score assumed.


Return: FASTQ file trimmed from the both ends (removed leading and trailing
		bases with quality lower than q)
'''

from Bio import SeqIO
import time


def get_data():
	file='/Users/Marlowe/Downloads/rosalind_bfil.txt'
	records=SeqIO.parse(file,"fastq")
# 	sequence_dict=SeqIO.to_dict(SeqIO.parse(file,"fastq"))

	with open(file) as file:
		quality=int(float(file.readline()))
		
	return quality,records

def count_from_middle(quality,records):
	quality_scores=[]
	sub_records=[]
	for record in records:
	# 	print(len(record.letter_annotations["phred_quality"]))
		'''count from the middle of the record down to find the first
			base with a q lower than quality'''
		start_base=0
		stop_base=0
		for i in range(len(record.letter_annotations["phred_quality"])//2,0,-1):
			if record.letter_annotations["phred_quality"][i]<quality:
				start_base=i
				break
		'''count from the middle of the record up to find the first
			base with a q lower than quality'''
		for t in range(len(record.letter_annotations["phred_quality"])//2,
						len(record.letter_annotations["phred_quality"]),1):
	# 		print('counting up to %i, at %i' % (len(record.letter_annotations["phred_quality"]),t))
	# 		print('testing
			if record.letter_annotations["phred_quality"][t]<quality:
	# 			print('found stop base')
				stop_base=t
				break
		'''create a sub_record of the current record starting and stopping at locations found above
		'''
		sub_records.append(record[start_base:stop_base])
	# 	print(sub_record)
	# 	print(record.letter_annotations["phred_quality"])
	# 	print(quality)
	# 	print('sub record start base: %i' % start_base)
	# 	print('sub record stop base: %i' % stop_base)
	# 	print('len of record: %i' % len(record.seq))
	# 	print('len of sub record: %i' % len(sub_record.seq))
	# 	print(record.letter_annotations["phred_quality"][0])
	# 	quality_scores.append(record.letter_annotations["phred_quality"])
	return sub_records

def remove_just_first(quality,records):
	quality_scores=[]
	sub_records=[]
	for record in records:
		start_base=0
		stop_base=0
# 		print('quality scores for current record= '+str(record.letter_annotations["phred_quality"]))
		'''count from the beginning of the record up to find the first
			base with a q higher than quality'''
		for i in range(0,len(record.letter_annotations["phred_quality"])):
			if record.letter_annotations["phred_quality"][i]>quality:
				start_base=i
# 				print('start_base: %i' % start_base)
				break
		'''count from the end of the record down to find the first
			base with a q higher than quality'''
		for t in range(len(record.letter_annotations["phred_quality"])-1,0,-1):
# 			print('counting down from %i, at %i' % (len(record.letter_annotations["phred_quality"]),t))
			if record.letter_annotations["phred_quality"][t]>quality:
				stop_base=t
# 				print('stop_base: %i of %i total bases' % (stop_base,len(record.letter_annotations["phred_quality"])-1))
				break
		'''create a sub_record of the current record starting and stopping at locations found above
		'''
# 		print('ammended quality scores: '+str(record[start_base:stop_base+1].letter_annotations["phred_quality"]))
		sub_records.append(record[start_base:stop_base+1])
	return sub_records			
				
def save_data(records):
	with open('text.txt','w') as newfile:
		for record in records:
			newfile.write(record.format("fastq"))
# 		for i in range(0,len(output_list)):
# 			writer.writerow(output_list[i])
		newfile.close()

def main():
	start_time = time.time()
	q,data=get_data()
# 	sub_records=count_from_middle(q,data)
	sub_records=remove_just_first(q,data)
	save_data(sub_records)
# 	for record in sub_records:
# 		print(record.format("fastq"))
# 	print(q)
	print('%f seconds' % (time.time() - start_time))
	
main()