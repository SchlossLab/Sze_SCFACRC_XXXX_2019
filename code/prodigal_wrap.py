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


############################################################################################

# Function to split files into small files if needed
def split_files(fasta_file, file_size_threshold, split_lines_threshold, temp_folder_name):
	# Split file if it is larger than file_size_threshold (default = 25MB)
	if os.path.getsize(fasta_file) > file_size_threshold:

		print("Input fasta is larger than %s bytes. Splitting..." % (file_size_threshold))

		os.system("split --suffix-length=7 --lines=%s %s %s/tmp_prodigal_" % 
			(split_lines_threshold, fasta_file, temp_folder_name))

	else:

		print("File is good to go and does not need splitting")

		os.system("cp %s %s" % (fasta_file, temp_folder_name))



# Function to run the prodigal call
def run_prodigal(temp_folder_name):

	if os.path.exists("%s/prodigal_run" % (workdir)):

		print("prodigal_run directory exists")

	else:

		print("No prodigal_run directory found, creating...")

		os.system("mkdir %s/prodigal_run" % (workdir))

	# Run prodigal call with standard defaults
	os.system("ls %s/* | xargs -I {} --max-procs=8 prodigal -q -c -i {} -o {}.genes \
-a {}.out -d {}.nucl -p meta" % (temp_folder_name))

	# Combine all files together and move to prodigal run directory

	os.system("%s/*.out > %s/prodigal_run/tmp_genes.fasta" % (temp_folder_name, workdir))

	os.system("%s/*.nucl > %s/prodigal_run/tmp_genes.nucleotide" % 
		(temp_folder_name, workdir))





# Run the actual program
def main(contig_fasta, FileSize, splitLines, temp_folder):
	
	split_files(contig_fasta, FileSize, splitLines, temp_folder)

	run_prodigal(temp_folder)
	
	
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
		Specify with the '-td' tag with the defaul set to 'tmp'.\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.contig_file, args.max_size, args.split_size, args.temp_dir)