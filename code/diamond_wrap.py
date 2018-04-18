#!python

# This is re-written python code to accomplish: 
	# Generate ORFs counts 

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

# Create a working directory for the analysis
def make_diamond_directory():

	# Create a directory to hold diamond analysis
	if os.path.exists("%sdiamond_analysis" % (workdir)):

		print("diamond_analysis directory exists")

	else:

		print("No diamond_analysis directory found, creating...")

		os.system("mkdir %sdiamond_analysis" % (workdir))



# Function to modify the fasta file
def make_new_fasta(opf_fasta_name):

	# read in mmseqs2 fasta file
	temp_file = open(opf_fasta_name, 'r')

	sequence_dict = {}
	sequence_name = []
	protein_seq = []
	
	print("Reading in fasta file %s" % (opf_fasta_name))
	
	# Create fasta dictionary
	for line in temp_file:

		if ">" in line and "start_type" in line:

			sequence_name.append(line)

		if ">" not in line and "start_type" not in line:

			protein_seq.append(line)

	# Create fasta dictionary
	for i, j in enumerate(sequence_name):

		sequence_dict[j] = protein_seq[i]

	# Create a new fasta file
	new_fasta = open("%sdiamond_analysis/temp_ref.fasta" % (workdir), 'w')

	print("Writing out new temp fasta file for diamond.")

	for protein_seq_name in sequence_dict:

		new_fasta.write(protein_seq_name+sequence_dict[protein_seq_name])

	new_fasta.close()


# Create a diamond based database
def make_diamond_db():
	
	os.system("diamond makedb --in %sdiamond_analysis/temp_ref.fasta \
-d %sdiamond_analysis/reference_db" % 
	 	(workdir, workdir))


# Run diamond of each fastq against the created database 
def run_diamond_alignment(sample_list):

	for i in sample_list:

		print("Running diamond analysis on sample %s" % (i))


		# Use blast to get hits of ORFs to Uniprot genes
		os.system("diamond blastx -q %s%s_hrm_r2.fastq -d %sdiamond_analysis/reference_db.dmnd \
-a %sdiamond_analysis/%s_output.daa -t %sdiamond_analysis/ \
--max-target-seqs 1 --evalue 1e-15 --id 0.90 --threads 8" % 
			(workdir, i, workdir, workdir, i, workdir))

		os.system("diamond view -a %sdiamond_analysis/%s_output.daa \
-o %sdiamond_analysis/%s.diamondout" % 
			(workdir, i, workdir, i))

		# Diamond output file has the following default columns (left to right)
		# qseqid = Query Seq - id, sseqid = Subject Seq - id, 
		# pident = Percentage of identical matches,  length = Alignment length, 
		# mismatch = Number of mismatches, gapopen = Number of gap openings, 
		# qstart = Start of alignment in query, qend = End of alignment in query, 
		# sstart = Start of alignment in subject, send = End of alignment in subject, 
		# evalue = Expect value, bitscore = Bit score

# Generate counts
def count_diamondoutput(sample_list):

	for i in sample_list:

		print("Obtaining orf counts from %s." % (i))

		os.system("cut -f 2 %sdiamond_analysis/%s.diamondout | sort | uniq -c | \
sed 's/^ *//' | sed 's/ /\\t/' | sed 's/$/\\t%s/' >> %sdiamond_analysis/orf_abund.tsv" % 
			(workdir, i, i, workdir))





# Run the actual program
def main(opf_fasta, sample_names):
	
	make_diamond_directory()

	make_new_fasta(opf_fasta)

	make_diamond_db()

	samples_to_be_used = create_samples_to_download(sample_names)

	run_diamond_alignment(samples_to_be_used)

	count_diamondoutput(samples_to_be_used)



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
	parser.add_argument("-s", "--sample_list", 
		default="%swhole_metagenome_samples.txt" % (workdir), 
		type=str, help="Text file with list of samples.\n")

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.opf_cluster_fasta, args.sample_list)















