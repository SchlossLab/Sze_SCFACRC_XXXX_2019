#!python

# This is re-written python code to accomplish two tasks
	# Generate ORFs
	# Get OPFs
# Translated Perl from Hannigan ProdigalWrapperLargeFiles.sh and ClusterOPFs.sh program

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Set up working directory
workdir = "data/raw/"

# Set up refrence directory
refdir = "data/references/"


############################################################################################

# Create a diamond based database
def make_diamond_db(opf_fasta_name):
	# Create a directory to hold diamond analysis
	if os.path.exists("%sdiamond_analysis" % (workdir)):

		print("diamond_analysis directory exists")

	else:

		print("No diamond_analysis directory found, creating...")

		os.system("mkdir %sdiamond_analysis" % (workdir))



	# read in mmseqs2 fasta file
	temp_file = open("test_head.fasta", 'r')

	sequence_dict = {}
	sequence_name = []
	protein_seq = []
	
	
	for line in temp_file:

		if ">" in line and "start_type" in line:

			sequence_name.append(line)

		if ">" not in line and "start_type" not in line:

			protein_seq.append(line)


	for i, j in enumerate(sequence_name):

		sequence_dict[j] = protein_seq[i]

	
		


	# os.system("diamond makedb --in %s -d %sdiamond_analysis/reference_db" % 
	# 	(opf_fasta_name, workdir))








# Run the actual program
def main(opf_fasta):
	
	make_diamond_db(opf_fasta)



# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-ocf", "--opf_cluster_fasta", 
		default="%smmseq2_opf_run/clu_seq.fasta" % (workdir), 
		type=str, help="Fasta file after mmseqs2 clustering. \
		This can be called with 'ocf' tag and the default is \
		%smmseq2_opf_run/clu_seq.fasta.\n" %(workdir))

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.opf_cluster_fasta)















