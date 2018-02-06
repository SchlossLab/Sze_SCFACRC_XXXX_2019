#!python

# This code is used to create a master file with all contigs.
# The step itself can be done using any language, I have chosen python to make
# it in part more easily readable.


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################

# Function to rename the given sequence file
def get_fasta_header(sampleList, contig_fa):

	# Set up empty dictionary
	temp_dict = {}


	for srr in sampleList:
		# set up read in by line
		fasta_file = open("%s%s/%s" % (workdir, srr, contig_fa), 'r')

		# set up empty variables and dictionaries to be used later
		contig_name = []
		contig = []
		
		# Read in each line of the stored fasta file
		for line in fasta_file:
			
			if line.find(">", 0) == 0:
				# removes the > and strips extraneous components
				temp_name = re.sub('>', ">%s_" % (srr), line)
				contig_name.append(temp_name.strip('\t\n'))
		
			# If it doesn't have the ">" assumes it is a sequence
			else: 
				# strips new line
				contig.append(line.strip('\n'))

		# Close the reading of the file
		fasta_file.close()

		# Create a dictionary with the sequence name and corresponding sequence
		for i, temp_name in enumerate(contig_name):
			# creates the actual dictionary of seq names and sequences
			temp_dict[temp_name] = contig[i]

	return(temp_dict)


#Function to write out the the combined dictionary into an overall fasta file
def create_complete_fasta(combinedDict):
	# Set up initial file call and directory variable to search
	file_dirs = os.listdir(workdir)
	file_present = False
	# Set initial counter to 0
	x = 0
	# iterate through the directory list and only stop if the file end is reached
	# or if file with the same name is found
	while file_present == False and x <= (len(file_dirs)-1):

		# check for file in directory
		if "all_contigs.fasta" in file_dirs[x]:
			# change to true if present
			file_present = True

		else:

			file_present = False
		# add one to the counter
		x += 1
	# Only create combined contig fasta if file doesn't already exist
	if file_present == False:
		# Sets up the file that will be read into
		test = open("%sall_contigs.fasta" % (workdir), 'w')

		for seq_name in combinedDict:
			temp_seq = combinedDict[seq_name]

			test.write(seq_name+'\n'+temp_seq+'\n')

		test.close() 

	else:
		print("File already exists.") 


# Runs the overall program 
def main(sampleListFile, contigNameEnding):

	samples_to_be_used = create_samples_to_download(sampleListFile)
	all_contig_dict = get_fasta_header(samples_to_be_used, contigNameEnding)
	create_complete_fasta(all_contig_dict)


# Initializes at the start of the program
if __name__ == '__main__': 

	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-s", "--sample_list", 
		default="%swhole_metagenome_samples.txt" % (workdir), 
		type=str, help="Text file with list of samples to run through the program.\n")
	parser.add_argument("-e", "--file_ending", 
		default="final.contigs.fa", type=str, 
		help="Ending assigned for each individual contig file. \
		The default is 'final.contigs.fa'.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.file_ending)