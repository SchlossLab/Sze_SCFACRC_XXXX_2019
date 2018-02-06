#!python

# This is re-written python code to accomplish two tasks
	# Convert SRA to fastq
	# Run quality control over each read per sample

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Set up working directory
workdir = "data/raw/"

# Set up shortcuts for keys and other download components
# These are changed to match Geoff's SRA data
asp = "/home/marcsze/.aspera/cli/etc/asperaweb_id_dsa.openssh"
sra_nih = "anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra"
database = "SRP"
seq_dir = "SRP108"
study = "SRP108915"

############################################################################################

# Function that takes command line arguments
def command_line():
	commands = sys.argv
	metagenome_sample_file = commands[1]

	return metagenome_sample_file

# Function to read in samples in stated text file for samples of interest
def create_samples_to_download(keep_file_name):

	temp_samples = []
	temp_file = open(keep_file_name, 'r')

	for line in temp_file:

		temp_samples.append(line.strip('\n'))

	# Close the reading of the file
	temp_file.close()	

	return(temp_samples)


# Cycles through and downloads each sample in d_files and moves it to appropriate directory
def download_files(keep_files):

	for files in keep_files:

		# Runs the actual download
		os.system("ascp -QTr -k 1 -l 1G -i %s %s/%s/%s/%s/%s/ %s" % 
			(asp, sra_nih, database, seq_dir, study, files, workdir))
		# Moves file back up a directory
		os.system("mv %s%s/%s.sra %s" % (workdir, files, files, workdir))
	# Removes extra empty directory
		os.system("rm -r %s%s" % (workdir, files))
	
		print("Completed downloading and moving sample: %s" % (files))



# Convert the sra file to fastq files (one forward (F) and one reverse (R))
def convert_sra_to_fastq(keep_files):

	for files in keep_files:

		print("Started creating fastq for: %s.sra" % (files))

		os.system("fastq-dump --split-3 %s%s.sra -O %s" % (workdir, files, workdir))

		print("Completed creating fastq for: %s.sra" % (files))


### Need to modify the cutadapt step based on jenior blog if seqs not from
### the SRA. 

# Function to run the cutadapt step 
def run_cutadpat(keep_files, metafileName, runReverseComplement):

	# set a counter
	x = 0
	# Set up the storage dictionaires
	F_three_temp_samples = {}
	R_three_temp_samples = {}
	F_five_temp_samples = {}
	R_five_temp_samples = {}
	# Open the needed data file (manually set above)
	temp_file = open(metafileName, 'r')
	# run through each line adding the forward and reverse seq components
	# to the respective dictionary
	for line in temp_file:
		# skip the header line that is read in
		if x >= 1:
			# remove the new line character
			temp_line = line.strip('\n')
			# split into a list based on the \t 
			temp_vector = temp_line.split('\t')
			# add a full sequence primer that includes adapter and index
			F_five_temp_samples[temp_vector[0]] = temp_vector[3]

			# temporarily store the reverse sequence primer
			F_three_temp_samples[temp_vector[0]] = temp_vector[4]
			
			temp_R_seq = temp_vector[4]
			temp_F_seq = temp_vector[3]

			if runReverseComplement == False:

				R_five_temp_samples[temp_vector[0]] = temp_vector[5]

				R_three_temp_samples[temp_vector[0]] = temp_vector[6]

			else:
				R_five_temp_samples[temp_vector[0]] = run_rev_complement(temp_R_seq)

				R_three_temp_samples[temp_vector[0]] = run_rev_complement(temp_F_seq)
		
				
		# Move counter up by 1
		x += 1
	# for each sample run the cutadapt code to remove them
	for sample_fastq in keep_files:
		# Runs the actual cutadapt command
		os.system("cutadapt --error-rate=0.1 --overlap=10 \
-a %s -g %s -A %s -G %s -o %s%s_adp_trim_1.fastq -p %s%s_adp_trim_2.fastq \
%s%s_qf_1.fastq %s%s_qf_2.fastq" % 
			(F_five_temp_samples[sample_fastq], 
				R_five_temp_samples[sample_fastq], 
				F_three_temp_samples[sample_fastq], 
				R_three_temp_samples[sample_fastq], 
				workdir, sample_fastq, workdir, sample_fastq, 
				workdir, sample_fastq, workdir, sample_fastq))



# Function that runs a reverse complement of a DNA string
def run_rev_complement(dna_sequence):

	# create an empty character place holder
	rev_seq = ""
	# run through each chr and change to the complement base
	for i in dna_sequence:

		if i == "A":
			rev_seq = "%sT" % (rev_seq)
		elif i == "T":
			rev_seq = "%sA" % (rev_seq)
		elif i == "C":
			rev_seq = "%sG" % (rev_seq)
		elif i == "G":
			rev_seq = "%sC" % (rev_seq)
		else:
			rev_seq = "%sN" % (rev_seq)

	# add the reverse complement 
	temp_rev_sequence = rev_seq[::-1]

	return(temp_rev_sequence)
	


# Quality filter and toss out reads that don't make the cut 
# -t quality filter of nucleotide (trimmed from the end of sequence)
# -Q average quality cutoff of sequence
# -l minimum length of sequence if shorter after trimming it is tossed
def run_quality_filter(keep_files, runningCutadapt):

	if runningCutadapt == False:

		for files in keep_files:

			print("Started quality filtering of F(1) and R(2) sample %s fastq" % 
				(files))

			os.system("sickle pe -f %s%s_1.fastq -r %s%s_2.fastq \
-t sanger -o %s%s_qf_1.fastq -p %s%s_qf_2.fastq \
-s %s%s_qf.orphan.fastq -q 30 -l 75" % 
			(workdir, files, workdir, files, workdir, files, 
				workdir, files, workdir, files))

			print("Completed quality filtering of sample %s" % (files))

	else:

		for files in keep_files:

			print("Started quality filtering of F(1) and R(2) sample %s fastq" % 
				(files))

			os.system("sickle pe -f %s%s_adp_trim_1.fastq -r %s%s_adp_trim_2.fastq \
-t sanger -o %s%s_qf_1.fastq -p %s%s_qf_2.fastq \
-s %s%s_qf.orphan.fastq -q 30 -l 75" % 
			(workdir, files, workdir, files, workdir, files, 
				workdir, files, workdir, files))

			print("Completed quality filtering of sample %s" % (files))



# Run the actual program
def main(sampleListFile, needCutadapt, primerfileName, needRC):
	samples_to_download = create_samples_to_download(sampleListFile)
	download_files(samples_to_download)
	convert_sra_to_fastq(samples_to_download)
	
	if needCutadapt == True:

		run_cutadpat(samples_to_download, primerfileName, needRC)
		run_quality_filter(samples_to_download, needCutadapt)

	else:

		run_quality_filter(samples_to_download, needCutadapt)
	
	
# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-s", "--sample_list", 
		default="%swhole_metagenome_samples.txt" % (workdir), 
		type=str, help="Text file with list of samples.\n")
	parser.add_argument("-ca", "--cut_adapt", default=False, 
		type=bool, help="Whether samples need adapters cut or not. \
		Default is 'false'. If 'true' need to provide a meta data file \
		with the '-m' tag or else program will fail.\n")
	parser.add_argument("-m", "--meta_data", 
		default="data/process/sequence_meta_data.tsv", 
		type=str, help="Must have the following order for columns in a \
		tab-delimited text file: sample_name, subjectID, diseaseClass, \
		five_prime_adapter, three_prime_adapter. \
		OPTIONAL to have columns that are also \
		rc_five_prime_adapter, rc_three_prime_adapter.\n")
	parser.add_argument("-rc", "--reverse_complement", 
		default=False, 
		type=bool, help="Argument to select whether to reverse complement \
		the provided reverse adapter. The default 'False' assumes that the data \
		provided is already reverse complemented.\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.cut_adapt, args.meta_data, args.reverse_complement)