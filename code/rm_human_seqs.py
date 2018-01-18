#!python

# This code is used to remove human seqs from fastq files.
# It uses samtools to complete this task but any program that IDs human sequences can be used


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"


############################################################################################





def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)


	print(samples_to_be_used)




if __name__ == '__main__': main()