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

# Function to generate a list of samples to grab based on dunns comparisons
# The sample id that matches clu_seq.fasta must be the first column
def generate_dunn_list(dunns_path_to_file):
	# Set up temp storage to only take unique samples
	temp_samples = set()
	temp_file = open(dunns_path_to_file, 'r')

	# Set up counter
	x = 1
	# Parse the csv file and take only the names
	for line in temp_file:
		# do not take the header line of the file
		if x > 1:

			test = line.split(',')

			temp_samples.add(test[0])

		# Move counter up
		x += 1

	
	return(list(temp_samples))


# Function to generate a dictionary with the seq_id and protein sequence
def generate_protein_dict(path_to_protein_fasta):
	# Set up temporary storage files
	temp_names = []
	temp_seq = []
	temp_dict = {}

	# Read in the protein fasta file
	temp_file = open(path_to_protein_fasta, 'r')
	# Grab the sequence name
	for line in temp_file:

		if ">" in line:

			if len(line) > 30:

				temp_value = line.split(" ")[0]

				temp_names.append(temp_value.strip('>'))

		# Grab the actual sequence
		else:

			temp_seq.append(line.strip('\n'))

	# Create the dictionary
	for pos, seq_name in enumerate(temp_names):


		temp_dict[seq_name] = temp_seq[pos]


	return(temp_dict)



# Create actual picking function and create new dictionary
def pick_fasta_sequence(seq_of_int, protein_db):
	# Set up storage variable
	temp_dict = {}
	# Go through the select list and pull sequences and create new dictionary
	for sample_name in seq_of_int:

		temp_sample = protein_db[sample_name]

		temp_dict[sample_name] = temp_sample


	return(temp_dict)


# Create a new pared down protein fasta file and write it out
def create_new_fasta(pared_down_dict, write_file_path):

	test = open(write_file_path, 'w')

	for seq_name in pared_down_dict:
			
		temp_seq = pared_down_dict[seq_name]

		test.write('>'+seq_name+'\n'+temp_seq+'\n')

	test.close()

	print("Completed pulling and creating protein sequence fasta.")



# Run the actual program
def main(dunnsFile, ProteinFastaFile, outputFilePath):
	
	dunn_list = generate_dunn_list(dunnsFile)

	protein_dict = generate_protein_dict(ProteinFastaFile)

	pared_protein_dict = pick_fasta_sequence(dunn_list, protein_dict)

	create_new_fasta(pared_protein_dict, outputFilePath)



	
	
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
	parser.add_argument("-of", "--output_file", 
		default="data/process/sig_dunn_protein.fasta", 
		type=str, help="This is the name of the outputed \
		protein fasta file based on the previous dunns test. \
		The default is 'data/process/sig_dunn_protein.fasta'. \n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.dunns_file, args.protein_fasta_file, args.output_file)