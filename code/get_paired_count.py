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


# Function to read in a fastq and count



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
	


if __name__ == '__main__': main()