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











# Run the actual program
def main():
	



# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-cf", "--contig_file", 
		default="%sall_contigs.fasta" % (workdir), 
		type=str, help="Fasta file with the assembled contigs.\n")
	parser.add_argument("-mfs", "--max_size", default=25000000, 
		type=int, help="Determines whether the contig file needs to be split \
		into smaller sub files for the program to run. \
		Specify with the '-mfs' tag with the default set at 25MB (25000000).\n")
	parser.add_argument("-ss", "--split_size", default=20, 
		type=int, help="Determines the total number of lines allowed in split file. \
		Specify with the '-ss' tag with the default set at '20'.\n")
	parser.add_argument("-td", "--temp_dir", default="%stmp" % (workdir), 
		type=str, help="Temporary directory in which split files will go. \
		Prodigal analysis will be completed within this temp dir. \
		Specify with the '-td' tag with the default set to 'tmp'.\n")
	parser.add_argument("-mmd", "--mmseqs2_dir", default="/home/marcsze/bin/MMseqs2", 
		type=str, help="Directoy on system the mmseqs2 program is stored in. \
		This will be used to set the path in linux based systems for the program call. \
		Specify with the '-mmd' tag with the default set to '/home/marcsze/bin/MMseqs2'.\n")
	parser.add_argument("-qd", "--query_database", 
		default="%sgenes.pep.format.fasta" % (refdir), 
		type=str, help="Reference data base used during OPF mmseq2 search. \
		This is used to specify the query fasta to make the query DB from. \
		Specify with the '-qd' tag with the default set to 'genes.pep.format.fasta' \
		from the KEGG 2016 genes database.\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main()















