#!python

# This code is used to generate a contig length table

# Translated Perl from Hannigan and should be run after the prodigal wrapper 
	# getOrfAbundance.sh program

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"

# Set up refrence directory
refdir = "data/references/"


############################################################################################

# Get the lengths that are needed
def get_contig_lengths(contig_fasta_name):
	# Read in contig fasta file
	fasta_file = open(contig_fasta_name, 'r')
	# inialize an empty dictionary
	contig_length_dict = {}
	# create the contig length dictionary
	print("Reading in contig fasta file: %s" % (contig_fasta_name))

	for line in fasta_file:

		if ">" in line:

			temp_lengths = line.split()[3]

			temp_names = line.split()[0]

			contig_length_dict[temp_names.strip(">")] = temp_lengths.strip("len=")

	
	return(contig_length_dict)

			
# Create the contig length table
def create_contig_length_table(length_dictionary):

	final_table = open("%scontig_length_table.tsv" % (workdir), 'w')

	print("Creating contig length table.")


	for seq_name in length_dictionary:
		

		final_table.write(seq_name+'\t'+length_dictionary[seq_name]+'\n')

	final_table.close()






# Run the actual program
def main(contig_fasta_file):

	length_dict = get_contig_lengths(contig_fasta_file)

	create_contig_length_table(length_dict)
	




# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-cf", "--contig_fasta", 
		default="%sall_contigs.fasta" % (workdir), 
		type=str, help="Fasta file with all the contigs. \
		This can be called with 'cf' tag and the default is \
		%sall_contigs.fasta.\n" %(workdir))

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.contig_fasta)