#!python

# This code is used to generate the total sequences within a given fasta file.
# It assumes that the programs that have been used maintain sequence
# pairing when processing the R1 (forward) and R2 (reverse) fastq files


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################


# Function to read in a fastq and count unique sequences
# Default would be to collect total and then 
# Sequences after quality filtering steps
def get_seq_counts(sampleList, seq_ending):

	temp_dict = {}

	for fastq_name in sampleList:

		print("Counting sequences in %s%s" % (fastq_name, seq_ending))
		
		x = 0

		temp_file_qf = open("%s%s%s" % 
			(workdir, fastq_name, seq_ending), 'r')

		for line in temp_file_qf:

			if "@SRR" in line:

				x += 1

		temp_dict[fastq_name] = x

	return(temp_dict)





# Example code for file read in
#temp_samples = []
#	temp_file = open(keep_file_name, 'r')

#	for line in temp_file:

#		temp_samples.append(line.strip('\n'))

	# Close the reading of the file
#	temp_file.close()	

#	return(temp_samples)


# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	test_full = get_seq_counts(samples_to_be_used, "_1.fastq")
	test_qf = get_seq_counts(samples_to_be_used, "_qf_1.fastq")
	test_hrm = get_seq_counts(samples_to_be_used, "_hrm_r1.fastq")

	print(test_full)
	print(test_qf)
	print(test_hrm)

if __name__ == '__main__': main()