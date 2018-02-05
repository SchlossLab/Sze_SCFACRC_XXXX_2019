#!python

# This code is used on aligned sam files to generate the necessary
# files used with concoct. Uses samtools, bedtools, and picard.


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

# Function to create bam files from sam files
def create_bam_files(sampleList, contigName, outputName):

	# create fai file
	print("Creating fai file.")
	os.system("samtools faidx %s.fasta" % (contigName))
	print("Completed creation of fai file.")
	
	# Iterate through each individual sam file and convert to bam
	for sam_file in sampleList:
		# convert to bam 
		print("Converting %s to bam format." % (sam_file))

		os.system("samtools view -bt %s.fasta.fai %s%s_%s_bowtie.sam > %s%s_%s.bam" % 
			(contigName, workdir, sam_file, outputName, workdir, sam_file, outputName))

		# sort the bam file
		print("Sorting created bam file.")

		os.system("samtools sort -o %s%s_%s_sort.bam %s%s_%s.bam" % 
			(workdir, sam_file, outputName, workdir, sam_file, outputName))

		#index the bam file
		print("Indexing the sorted and created bam file.")
		
		os.system("samtools index %s%s_%s_sort.bam" % 
			(workdir, sam_file, outputName))

		print("Completed %s bam creation." % (sam_file))


# Function to call picard to mark and remove duplicates
def mark_duplicates(sampleList, outputName):

	for bam_file in sampleList:
		# mark and remove duplicates
		print("Marking and removing duplicate sequences from %s." % 
			(bam_file))

		os.system("java -jar $PICARD_JARS/picard.jar MarkDuplicates \
INPUT=%s%s_%s_sort.bam \
OUTPUT=%s%s_%s_sort_rmdup.bam \
METRICS_FILE=%s%s_%s_sort_rmdup.metrics \
AS=TRUE \
VALIDATION_STRINGENCY=LENIENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
REMOVE_DUPLICATES=TRUE" % 
(workdir, bam_file, outputName, workdir, bam_file, outputName, 
	workdir, bam_file, outputName))

		# Re sort the resulting bam file
		print("Re-sorting the bam file.")

		os.system("samtools sort -o %s%s_%s_sort_rmdup_sort.bam %s%s_%s_sort_rmdup.bam" % 
			(workdir, bam_file, outputName, workdir, bam_file, outputName))
		
		# Re index the resulting re-sorted bam file
		print("Re-indexing the re-sorted bam file.")

		os.system("samtools index %s%s_%s_sort_rmdup_sort.bam" % 
			(workdir, bam_file, outputName))


# Function to create final needed coverage and linkage files for concoct
def get_cover_and_link(sampleList, contigName, outputName):

	# execute coverage data command for every sample
#	for bam_file in sampleList:
		# get coverage data
#		print("Obtaining coverage data for %s." % (bam_file))

#		os.system("genomeCoverageBed -ibam %s%s_%s_sort_rmdup_sort.bam \
#-g %s.fasta > %s%s_%s_final_contigs_coverage.txt" % 
#(workdir, bam_file, outputName, contigName, workdir, bam_file, outputName))

	# create a temp text file with the sample names
	cov_file = open("%stemp_coverage_file_names.txt" % (workdir),'w')

	# add each sample to the temp text file
	print("Creating a temporary sample file for coverage analysis.")

	for sampleN in sampleList:

		cov_file.write("%s%s_%s_final_contigs_coverage.txt" % 
			(workdir, sampleN, outputName)+'\n')
	# close the file
	cov_file.close()
	
	# creates the coverage table
	print("Creating coverage table.")

	os.system("python2 /sw/med/centos7/concoct/0.4.1/bin/gen_input_table.py \
--isbedfiles \
--samplenames %stemp_coverage_file_names.txt \
%s.fasta \
%s*_final_contigs_coverage.txt > \
%soverall_coverage_table.tsv" % 
(workdir, contigName, workdir, workdir))


	# create a temp text file with the sample names
	link_file = open("%stemp_linkage_file_names.txt" % (workdir),'w')

	# add each sample to the temp text file
	print("Creating a temporary sample file for linkage analysis.")

	for bam_file in sampleList:

		link_file.write("%s%s_%s_sort_rmdup_sort.bam" % 
			(workdir, bam_file, outputName)+'\n')
	# close the file
	link_file.close()


	# creates the linkage table
	print("Creating linkage table.")

	os.system("python2 /sw/med/centos7/concoct/0.4.1/bin/bam_to_linkage.py \
-m 8 --regionlength 500 \
--fullsearch --samplenames %stemp_linkage_file_names.txt \
%s.fasta \
%s*_%s_sort_rmdup_sort.bam > %soverall_linkage_table.tsv" % 
(workdir, contigName, workdir, outputName, workdir))

	# Remove the temporary sample text file
	os.system("rm %stemp_coverage_file_names.txt" % (workdir))
	os.system("rm %stemp_linkage_file_names.txt" % (workdir))



# Runs the overall program 
def main(sampleListFile, contigFile, outputEnding):
	# read in file list from -s call
	samples_to_be_used = create_samples_to_download(sampleListFile)
	#create_bam_files(samples_to_be_used, contigFile, outputEnding)
	#mark_duplicates(samples_to_be_used, outputEnding)
	get_cover_and_link(samples_to_be_used, contigFile, outputEnding)


# Upon program call executes these commands automatically
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-l", "--sample_list", default="%swhole_metagenome_samples.txt" % (workdir), type=str, help="Text file with list of samples\n")
	parser.add_argument("-c", "--contig_file", default="%sall_contigs_1kbto10kb" % (workdir), type=str, help="Combined cut contig fasta fil name\n")
	parser.add_argument("-s", "--sam_file_end", default="contig_1to10kb", type=str, help="Unique Ending of Samfiles\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.contig_file, args.sam_file_end)