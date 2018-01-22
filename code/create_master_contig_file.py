#!python

# This code is used to create a master file with all contigs.
# The step itself can be done using any language, I have chosen python to make
# it in part more easily readable.


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
contig_fa = "final.contigs.fa"


############################################################################################

# Function to rename the given sequence file
def get_fasta_header(sampleList):

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
	
	test = open("%sall_contigs.fasta" % (workdir), 'w')

	for seq_name in combinedDict:
		temp_seq = combinedDict[seq_name]

		test.write(seq_name+'\n'+temp_seq+'\n')

	test.close()  





# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	all_contig_dict = get_fasta_header(samples_to_be_used)
	create_complete_fasta(all_contig_dict)

	#print(test[">SRR5665186_k101_48 flag=1 multi=6.0000 len=1359"])



if __name__ == '__main__': main()