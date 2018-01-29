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
contig_fa = "final.contigs.fa"


############################################################################################





# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	all_contig_dict = get_fasta_header(samples_to_be_used)
	create_complete_fasta(all_contig_dict)


if __name__ == '__main__': main()