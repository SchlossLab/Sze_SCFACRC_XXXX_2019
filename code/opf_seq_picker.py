#!python

# This code will take the dunn result csv file and pick specific
# protein sequences from a fasta file and ouput to a new file.

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Set up working directory
workdir = "data/raw/"

# Set up refrence directory
refdir = "data/references/"


############################################################################################


def generate_dunn_list(dunns_path_to_file):

	temp_samples = []
	temp_file = open(dunns_path_to_file, 'r')

	for line in temp_file:

		test = line.split(',')

		temp_samples.append(test[0])

	return(temp_samples)

















# Run the actual program
def main(dunnsFile, ProteinFastaFile):
	
	dunn_list = generate_dunn_list(dunnsFile)


	
	
# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-dcf", "--dunns_file", 
		default="data/process/opf_dunn_select_results.csv", 
		type=str, help="Dunn comparison file. \
		The default is 'data/process/opf_dunn_select_results.csv'. \n")
	parser.add_argument("-pff", "--protein_fasta_file", 
		default="%smmseq2_opf_run/clu_seq.fasta" % (workdir), 
		type=str, help="This is the protein fasta file to pull sequences from. \
		The default is 'data/raw/mmseq2_opf_run/clu_seq.fasta'. \n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.dunns_file, args.protein_fasta_file)