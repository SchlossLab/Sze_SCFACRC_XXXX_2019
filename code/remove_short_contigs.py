#! python

# This code is used to remove contigs from the master contig list that is 
# under 1kb long. 


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"

############################################################################################

# Function to read in the all_contig file and return need dictionaries
def get_needed_data(sampleList, all_contig):

	# Set up empty dictionaries
	temp_dict = {}
	count_dict = {}
	name_dict = {}

	contig_file = open(all_contig, 'r')

	for line in contig_file:
		
		if ">" in line:

			all_info = line.split(" ")

			temp_length = all_info[3].strip('len=\n')

			tempName = all_info[0]

			temp_full_name = line.strip('\n')

		else:

			sequence_info = line.strip('\n')

			temp_dict[tempName] = sequence_info

			count_dict[tempName] = temp_length

			name_dict[tempName] = temp_full_name


	return temp_dict, count_dict, name_dict


# Function to only write back out to file contigs larger than 1Kb
def remove_contigs(seqDict, countDict, nameDict, cutoff):

	tempName = {}
	tempSeq = {}

	for contig in nameDict:

		if int(countDict[contig]) > int(cutoff):

			tempName[contig] = nameDict[contig]

			tempSeq[contig] = seqDict[contig]


	return tempName, tempSeq

# Function to create a fasta with only 1Kb or above contigs
def create_new_contig_file(nameDict, contigDict, output_file):

	# Set up initial file call and directory variable to search
	file_dirs = os.listdir(workdir)
	file_present = False
	# Set initial counter to 0
	x = 0
	# iterate through the directory list and only stop if the file end is reached
	# or if file with the same name is found
	while file_present == False and x <= (len(file_dirs)-1):

		# check for file in directory
		if output_file in file_dirs[x]:
			# change to true if present
			file_present = True

		else:

			file_present = False
		# add one to the counter
		x += 1
	# Only create combined contig fasta if file doesn't already exist
	if file_present == False:
		# Sets up the file that will be read into
		test = open(output_file, 'w')

		for seq_name in contigDict:
			
			temp_name = nameDict[seq_name]
			temp_contig = contigDict[seq_name]

			test.write(temp_name+'\n'+temp_contig+'\n')

		test.close() 

	else:
		print("File already exists.") 



# Runs the overall program 
def main(sampleListFile, combinedContigFileName, min_contig_length, output_name):

	samples_to_be_used = create_samples_to_download(sampleListFile)
	seq_dict, length_dict, fullname_dict = get_needed_data(samples_to_be_used, combinedContigFileName)
	final_name_dict, final_contig_dict = remove_contigs(seq_dict, length_dict, fullname_dict, min_contig_length)
	create_new_contig_file(final_name_dict, final_contig_dict, output_name)
	


# Initializes at the start of the program
if __name__ == '__main__': 

	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-s", "--sample_list", 
		default="%swhole_metagenome_samples.txt" % (workdir), 
		type=str, help="Text file with list of samples to run through the program.\n")
	parser.add_argument("-cf", "--contig_file", 
		default="%sall_contigs.fasta" % (workdir), type=str, 
		help="The combined contig file to be parsed. \
		The default is 'all_contigs.fasta'.\n")
	parser.add_argument("-mc", "--min_contig", 
		default=1000, type=int, 
		help="The minimum size a contig needs to be to be kept. \
		The default is 1000 bp.\n")
	parser.add_argument("-o", "--output", 
		default="%sall_contigs_greater1kb.fasta" % (workdir), type=str, 
		help="File to save the remaining contigs to. \
		The default is 'all_contigs_greater1kb.fasta'.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.contig_file, args.min_contig, args.output)