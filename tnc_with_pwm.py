from Bio import SeqIO
import numpy as np
import re
import os 
import sys

triplet_to_base_array_index_mapping = dict()
max_seq_len = 0

def triplet_init_in_dictionary():
	counter = 0;
	base_values = ['A','T','C','G']

	for i in range(0,4):
		for j in range (0,4):
			for k in range (0,4):
				
				tnc = base_values[i]+ base_values[j]+ base_values[k]

				triplet_to_base_array_index_mapping[tnc] = counter

				counter+=1


def get_relative_frequency_base_sequence(base_sequence):

	frequency_based_count = np.zeros(64)

	triplet_count = len(base_sequence)-2;

	for i in range(triplet_count):
		tnc = base_sequence[i:i+3].upper()
		frequency_based_count[triplet_to_base_array_index_mapping[tnc]] += 1
	
	frequency_based_count = frequency_based_count/triplet_count
	return frequency_based_count

def write_feature_names(output_file_name):
	output_file_obj = open(output_file_name,"w")
	
	features = list(triplet_to_base_array_index_mapping.keys())
	for feature in features:
		output_file_obj.write(feature+',')
	output_file_obj.write('label\n')
	output_file_obj.close()


def create_position_probability_matrix(input_file_names):
	global max_seq_len

	temp_seq_list = []
	#counter = 0
	for input_file_name in input_file_names:		
		#print(input_file_name)
		for record in SeqIO.parse(input_file_name, "fasta"):
			sequence = re.sub('[0123456789]','',str(record.seq))
			temp_seq_list.append(sequence)
			if len(sequence)>max_seq_len:
				max_seq_len = len(sequence)
			
	position_probability_matrix  = np.zeros((64,max_seq_len))
	#print(counter)
	for sequence in temp_seq_list:
		for k in range(len(sequence)-2):
			tri_nucleotide = sequence[k:k+3].upper()
			position_probability_matrix[triplet_to_base_array_index_mapping[tri_nucleotide]][k] += 1

	position_probability_matrix = position_probability_matrix/len(temp_seq_list)

	return position_probability_matrix


def generate_ppm_based_feature(input_ppm, input_seq_mask):
	global max_seq_len
	feature_val_array = np.zeros(64)
	
	for k in range(64):
		diagonal = input_ppm[k,:] * input_seq_mask[:,k]
		feature_val_array[k] = np.sum(diagonal)

	return feature_val_array


def get_label(file_name, mode ):
	if mode == 0:
		if file_name.startswith('non'):
			return 0
		else:
			return 1

	else:
		if file_name.startswith('non'):
			return 0
		if file_name.startswith('strong'):
			return 1
		if file_name.startswith('weak'):
			return 2



def create_data_file(input_files,output_file,ppm):
	global max_seq_len
	write_feature_names(output_file)
	out = open(output_file,"a")

	for input_file in input_files:
		
		for record in SeqIO.parse(input_file, "fasta"):
			
			sequence = re.sub('[0123456789]','',str(record.seq))
			mask = np.zeros((max_seq_len,64))
			
			for k in range(len(sequence)-2):
				tnc = sequence[k:k+3].upper()
				mask[k][triplet_to_base_array_index_mapping[tnc]] = 1

			pr_data = generate_ppm_based_feature(ppm,mask)
			rf_data = get_relative_frequency_base_sequence(sequence)

			data = pr_data**2 + rf_data**2
			data = data**0.5
			if len(data) != 64:
				print('WTF')
			
			for k in range(64):
				out.write(str(data[k])+',')
			out.write(str(get_label(input_file,0))+'\n')

	out.close()

def create_data_file_for_relative_frequency(input_files,output_file):
	
	write_feature_names(output_file)
	out = open(output_file,"a")

	for input_file in input_files:
		
		for record in SeqIO.parse(input_file, "fasta"):
			
			sequence = re.sub('[0123456789]','',str(record.seq))
			
			data = get_relative_frequency_base_sequence(sequence)
			
			for k in range(64):
				out.write(str(data[k])+',')
			out.write(str(get_label(input_file,0))+'\n')

	out.close()

def main():
	triplet_init_in_dictionary()
	files = [sys.argv[1],sys.argv[2],sys.argv[3]]
	ppm = create_position_probability_matrix(files)
	print(ppm)
	create_data_file(files,sys.argv[4],ppm)
	#create_data_file_for_relative_frequency(files,sys.argv[4])
	print('Ready for training')
	

main()


