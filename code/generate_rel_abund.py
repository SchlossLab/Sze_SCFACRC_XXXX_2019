#!python

# This code is used to create a master bowtie reference file, aligning 
# sequences to contigs, and then using this reference generate 
# contig relative abundnace.
# This task was completed using BowTie2


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
all_contigs = "all_contigs.fasta"
reference = "bowtieReference"
summarydir = "data/process/tables/"

############################################################################################

# Function that creates a bowtie reference file from all contigs
def make_ref_file():

	print("Starting to build reference database")

	os.system("bowtie2-build -q %s%s %s/bowtieReference" % 
		(workdir, all_contigs, workdir))

	print("Completed building reference database")


# Function to align reverse sequences to the master contig reference database
def run_alignment(sampleList):

	for r2_fasta in sampleList:

		os.system("bowtie2 -x %s%s -q %s%s_qf_2.fastq -S %s%s_bowtie.sam -p 8 -L 25 -N 1" % 
			(workdir, reference, workdir, r2_fasta, workdir, r2_fasta))


# Function to create a dictionary with 0 counts for all the contigs 
def create_contig_dict():

	print("Generating names for contig dictionary")

	tempDict = {}

	temp_file = open("%s%s" % (workdir, all_contigs), 'r')

	for line in temp_file:

		if ">" in line:

			all_info = line.split(" ")

			temp_length = all_info[3].strip('len=\n')

			tempName = all_info[0].strip('>')

			tempDict[tempName] = temp_length


	temp_file.close()

	print("Completed generation of names list for contig dictionary")

	return tempDict



# Function to count specific abundances within the generated sam files
def get_contig_abundance(sampleList, contig_dict):

	countNames = contig_dict.keys()

	for sample_id in sampleList:

		tempDict = {}

		for sample_name in countNames:

			tempDict[sample_name] = 0


		temp_file = open("%s%s_bowtie.sam" % (workdir, sample_id), 'r')


		print("Generating count data for %s" %(sample_id))

		for line in temp_file:
			# regex line matching (taken from Geof perl code)
			if re.match('^@[A-Z]{2}', line) is None:

				all_info = line.split('\t')

				test = all_info[2]

				if test in tempDict:

					tempDict[test] += 1

		
		temp_file.close()	

		print("Writing count data for %s" %(sample_id))

		# Opens the necessary summary file
		write_file = open("%s%s_contig_rel_abund.tsv" % 
			(workdir, sample_id),'w')
		# Counter
		x = 0
		# Iterates through every sample 
		for sample_name in countNames:
			# Checks if it is the first sample
			if x == 0:
				# Adds header and sample counts
				write_file.write("sample_name"+'\t'+"contig"+'\t'+ \
					"contig_length"+'\t'+"count"+'\n'+ \
					sample_id+'\t'+sample_name+'\t'+ \
					str(contig_dict[sample_name])+'\t'+ \
					str(tempDict[sample_name])+'\n')

			else:
				# Adds only sample counts
				write_file.write(sample_id+'\t'+ \
					sample_name+'\t'+str(contig_dict[sample_name])+'\t'+ \
					str(tempDict[sample_name])+'\n')

			# adds to count
			x += 1
		# Clost the file
		write_file.close()

		print("Finished counting for %s" %(sample_id))
		

# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	#make_ref_file()
	#run_alignment(samples_to_be_used)	
	contig_name_dict = create_contig_dict()
	get_contig_abundance(samples_to_be_used, contig_name_dict)

if __name__ == '__main__': main()