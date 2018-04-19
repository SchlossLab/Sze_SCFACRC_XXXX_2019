#!python

# This is re-written python code to accomplish two tasks
	# Generate ORFs
	# Get the OPF groupings
# Translated Perl from Hannigan ProdigalWrapperLargeFiles.sh and ClusterOPFs.sh program

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Set up working directory
workdir = "data/raw/"

# Set up refrence directory
refdir = "data/references/"


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

	os.system("cat %s/*.out > %s/prodigal_run/tmp_genes.fasta" % (temp_folder_name, workdir))

	os.system("cat %s/*.nucl > %s/prodigal_run/tmp_genes.nucleotide" % 
		(temp_folder_name, workdir))


# Function to run mmseqs2 to establish the OPFs
def run_mmseqs2(mmseq2_dir_path, query_fasta):

	os.system("PATH=%s/bin:$PATH" % (mmseq2_dir_path))
	# Create mmseq2 directory
	if os.path.exists("%smmseq2_opf_run" % (workdir)):

		print("mmseq2_opf_run directory exists")

	else:

		print("No mmseq2_opf_run directory found, creating...")

		os.system("mkdir %smmseq2_opf_run" % (workdir))

	# Create a tmp directory within the mmseq2 directory
	if os.path.exists("%smmseq2_opf_run/tmp" % (workdir)):

		print("tmp directory in mmseq2_opf_run directory exists")

	else:

		print("No tmp directories found in mmseq2_opf_run, creating...")

		os.system("mkdir %smmseq2_opf_run/tmp" % (workdir))

		os.system("mkdir %smmseq2_opf_run/tmp_result" % (workdir))

	# Set up query DB
	os.system("mmseqs createdb %s %smmseq2_opf_run/query_DB" % (query_fasta, workdir))

	# Set up target DB
	os.system("mmseqs createdb %sprodigal_run/tmp_genes.fasta \
%smmseq2_opf_run/target_DB" % (workdir, workdir))

	# Create index 

	os.system("mmseqs createindex %smmseq2_opf_run/target_DB %smmseq2_opf_run/tmp" % 
		(workdir, workdir))

	# Create alignment

	os.system("mmseqs search %smmseq2_opf_run/query_DB %smmseq2_opf_run/target_DB \
%smmseq2_opf_run/resultDB %smmseq2_opf_run/tmp" % (workdir, workdir, workdir, workdir))


	# Create BLAST formatted file of Result DB (tsv based file)
	os.system("mmseqs convertalis %smmseq2_opf_run/query_DB %smmseq2_opf_run/target_DB \
%smmseq2_opf_run/resultDB %smmseq2_opf_run/resultDB.m8" % 
		(workdir, workdir, workdir, workdir))


	
	# Cluster the target DB
	os.system("mmseqs cluster %smmseq2_opf_run/target_DB %smmseq2_opf_run/target_clu \
%smmseq2_opf_run/tmp -e 0.001 --min-seq-id 0.4" % (workdir, workdir, workdir))


	# Create a target clustered sequence file
	os.system("mmseqs createseqfiledb %smmseq2_opf_run/target_DB \
%smmseq2_opf_run/target_clu %smmseq2_opf_run/clu_seq" % (workdir, workdir, workdir))


	# Create target clustered fasta file
	os.system("mmseqs result2flat %smmseq2_opf_run/target_DB %smmseq2_opf_run/target_DB \
%smmseq2_opf_run/clu_seq %smmseq2_opf_run/clu_seq.fasta" % 
		(workdir, workdir, workdir, workdir))


	# Create a clustered target sequence tsv file
	os.system("mmseqs createtsv %smmseq2_opf_run/target_DB %smmseq2_opf_run/target_DB \
%smmseq2_opf_run/target_clu %smmseq2_opf_run/clu.tsv" % 
		(workdir, workdir, workdir, workdir))








# Run the actual program
def main(contig_fasta, FileSize, splitLines, temp_folder, mmseq2_dir, query_DB_fasta):
	
	split_files(contig_fasta, FileSize, splitLines, temp_folder)

	run_prodigal(temp_folder)

	run_mmseqs2(mmseq2_dir, query_DB_fasta)
	
	
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
	main(args.contig_file, args.max_size, args.split_size, args.temp_dir, 
		args.mmseqs2_dir, args.query_database)