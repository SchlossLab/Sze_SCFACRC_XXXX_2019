#!python

# This code is used to create a master bowtie reference file, aligning 
# sequences to contigs, and then using this reference generate 
# contig relative abundnace.
# This task was completed using BowTie2


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
#all_contigs = "all_contigs.fasta"
#reference = "bowtieReference"
#summarydir = "data/process/tables/"

############################################################################################

# Function that creates a bowtie reference file from all contigs
def make_ref_file(contigPath, referencePath):

	print("Starting to build reference database")

	os.system("bowtie2-build -q %s %s" % 
		(contigPath, referencePath))

	print("Completed building reference database")


# Function to align reverse sequences to the master contig reference database
def run_alignment(sampleList, referencePath, outputEnding):

	for r2_fasta in sampleList:

		os.system("bowtie2 -x %s -q %s%s_qf_2.fastq -S %s%s_%s_bowtie.sam -p 8 -L 25 -N 1" % 
			(referencePath, workdir, r2_fasta, workdir, r2_fasta, outputEnding))


# Function to create a dictionary with the lengths of the contigs 
def create_contig_dict(contigPath, size_limit_used):

	print("Generating names for contig dictionary")

	tempDict = {}
	tempName = []
	temp_length = []

	temp_file = open("%s" % (contigPath), 'r')
	# Parse each line of the contig file
	for line in temp_file:
		# This is performed if using the defaults without any size selection
		if ">" in line and size_limit_used == False:

			all_info = line.split(" ")

			temp_length = all_info[3].strip('len=\n')

			tempName = all_info[0].strip('>')

			tempDict[tempName] = temp_length

		# From here beyond is done if the contig file has a size selection included
		elif ">" in line and size_limit_used == True:

			tempName.append(line.strip('>\n')) 

		elif ">" not in line and size_limit_used == True:

			x = 1

			for nucleotide in line:

				x +=1

			temp_length.append(x) 

	
	if size_limit_used == True:

		for i, contig in enumerate(tempName):

			tempDict[contig] = temp_length[i]

	print("Completed grabbing length of contigs and dictionary creation")

	temp_file.close()


	return(tempDict)



# Function to count specific abundances within the generated sam files
def get_contig_abundance(sampleList, contig_dict, outputEnd):

	countNames = contig_dict.keys()

	for sample_id in sampleList:

		tempDict = {}

		for sample_name in countNames:

			tempDict[sample_name] = 0


		temp_file = open("%s%s_all_contig_bowtie.sam" % (workdir, sample_id), 'r')


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
		write_file = open("%s%s_%s_rel_abund.tsv" % 
			(workdir, sample_id, outputEnd),'w')
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
def main(sampleListFile, contigFile, refFile, contigSizeLimit, outputEnding):
	# read in file list from -s call
	samples_to_be_used = create_samples_to_download(sampleListFile)
	# # create a bowtie reference with the -c and -r input calls
	make_ref_file(contigFile, refFile)
	# # # run alignment call with -r and -o input calls
	run_alignment(samples_to_be_used, refFile, outputEnding)	
	# generate contig names based on the -c input call
	contig_length_dict = create_contig_dict(contigFile, contigSizeLimit)
	# create contig abundance tables with the -o input call
	get_contig_abundance(samples_to_be_used, contig_length_dict, outputEnding)

# Upon program call executes these commands automatically
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-s", "--sample_list", default="%swhole_metagenome_samples.txt" % (workdir), type=str, help="Text file with list of samples\n")
	parser.add_argument("-c", "--contig_file", default="%sall_contigs.fasta" % (workdir), type=str, help="Combined contig fasta file\n")
	parser.add_argument("-r", "--reference", default="%sbowtieReference" % (workdir), type=str, help="Bowtie2 Reference file name\n")
	parser.add_argument("-sl", "--size_limit", default=False, 
		type=bool, help="Whether contig fasta file is size limited or not. \
		The default is that there is no size limitation (e.g. 'False'). \n")
	parser.add_argument("-o", "--output", default="all_contig", type=str, help="Universal output file ending (outputs in tab format)\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.contig_file, args.reference, args.size_limit, args.output)