#!python

# This code selectively filters for specific KEGG IDS 
# Designed to search the orf_gene_alignment.tsv file for OPFs that match the given KEGG IDs
# This was used to pull out only the specific genes for butyrate, propionate, and acetate


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

# Function to read in the text file with KEGG IDs of interest
def get_ids_of_int(keggids):

	print("Reading in the specific KEGG IDs.")

	temp_file = open(keggids, 'r')

	temp_list = []

	for kegg in temp_file:

		temp_list.append(kegg.strip('\n'))

	return(temp_list)


# Function to create a dictionary of the kegg to OPF calls
def make_opf_dict(orf_path):

	print("Creating the search dictionary")

	temp_file = open(orf_path, 'r')

	temp_contig_dict = {}

	for opf in temp_file:

		temp_line = opf.split('\t')

		temp_kegg_list = temp_line[0].split("|")

		temp_contig_id = temp_line[1]

		temp_kegg_id = temp_kegg_list[len(temp_kegg_list)-1]

		temp_contig_dict[temp_contig_id] = temp_kegg_id

	return(temp_contig_dict)



# Function to perform the match throughout the dictionary and keep those that match
def make_matches(list_of_keggs, dict_of_contigs):

	print("Finding matches from specified file.")

	temp_dict = {}

	for contig in dict_of_contigs:

		temp_kegg = dict_of_contigs[contig]

		if temp_kegg in list_of_keggs:

			temp_dict[contig] = temp_kegg


	return(temp_dict)


# Function to write out matched data into a new tsv file
def write_the_file(final_dict, output_name):

	write_file = open("%s" % (output_name),'w')

	for contig in final_dict:

		write_file.write(contig+'\t'+str(final_dict[contig])+'\n')


	print("Wrote match results to %s" % (output_name))




# Run the actual program
def main(keggids_path, orf_alignment_data, outputName):

	keggIDs = get_ids_of_int(keggids_path)

	contig_dict = make_opf_dict(orf_alignment_data)

	matched_dict = make_matches(keggIDs, contig_dict)

	write_the_file(matched_dict, outputName)




# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-iko", "--input_kegg_ids", 
		default="data/process/scfa_kegg_ids.txt", 
		type=str, help="This is a text file that contains the KEGG IDs \
		of interest. Each line should contain a unique ID. \
		The default file used is 'data/process/scfa_kegg_ids.txt'.\n")
	parser.add_argument("-gad", "--gene_alignment_data", 
		default="data/process/orf_gene_alignment.tsv", 
		type=str, help="This is a tab separated file that contains the KEGG \
		matches for each specific contig. \
		The default file is 'data/process/orf_gene_alignment.tsv'.\n")
	parser.add_argument("-o", "--output_file", 
		default="data/process/select_scfa_opf_matches.tsv", 
		type=str, help="A tsv file of all the OPFs that match to the KEGG IDs of interst. \
		The default file used for the output to be saved under is \
		'data/process/select_scfa_opf_matches.tsv'.\n")

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.input_kegg_ids, args.gene_alignment_data, args.output_file)
