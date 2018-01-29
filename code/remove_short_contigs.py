#! python

# This code is used to remove contigs from the master contig list that is 
# under 1kb long. 


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
all_contig = "all_contigs.fasta"
output_file = "all_contigs_greater1kb.fasta"


############################################################################################

# Function to read in the all_contig file and return need dictionaries
def get_needed_data(sampleList):

	# Set up empty dictionaries
	temp_dict = {}
	count_dict = {}

	contig_file = open("%s%s" % (workdir, all_contig), 'r')

	for line in contig_file:
		
		if ">" in line:

			all_info = line.split(" ")

			temp_length = all_info[3].strip('len=\n')

			tempName = all_info[0]

		else:

			sequence_info = line

			temp_dict[tempName] = sequence_info

			count_dict[tempName] = temp_length


	return temp_dict, count_dict


# Function to only write back out to file contigs larger than 1Kb



# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	seq_dict, length_dict = get_needed_data(samples_to_be_used)

	print(len(seq_dict), len(length_dict))



if __name__ == '__main__': main()