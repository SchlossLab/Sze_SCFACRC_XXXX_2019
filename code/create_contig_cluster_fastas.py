#!python

# This code is designed to take a contig fasta and CONCOCT cluster table and
# sort them into clustered fasta files for later searchs against 
# sequence databases of your choice.



############# Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################

# Function to create needed dictionaries
def create_dictionaires(contig_fasta_path, cluster_table_path):

	print("Creating needed dictionaries.")

	temp_contig_dict = {}
	temp_cluster_dict = {}

	temp_fasta = open(contig_fasta_path, 'r')

	for line in temp_fasta:

		if ">" in line:

			temp_name = line.strip('>\n')

		else:

			temp_contig_dict[temp_name] = line.strip('\n')


	temp_fasta.close()

	temp_cluster = open(cluster_table_path, 'r')

	for line in temp_cluster:

		temp_info = line.split(',')

		temp_name = temp_info[0]

		temp_group = temp_info[1].strip('\n')

		temp_cluster_dict[temp_name] = int(temp_group)


	temp_cluster.close()

	return(temp_contig_dict, temp_cluster_dict)


# Function to create a unique list of the total cluster groups
def get_cluster_groups(cluster_data):

	print("Getting cluster groups.")

	temp_clusters = []

	for seq in cluster_data:

		temp_group = cluster_data[seq]

		if temp_group not in temp_clusters:

			temp_clusters.append(temp_group)



	return(temp_clusters)




# Function to create the needed output directory to put fasta files
def create_output_dir(outputDir):

	# Create a directory to hold diamond analysis
	if os.path.exists(outputDir):

		print("Clustered fasta directory exists")

	else:

		print("No clustered fasta directory found, creating...")

		os.system("mkdir %s" % (outputDir))





# Function to selectively add sequences to specific fasta files
def create_clustered_fastas(cluster_info, contigDict, clusterDict, output):

	for i in cluster_info:

		print("Making fasta file for cluster group %s." % (i))

		test_file = open("%scluster_group_%s.fasta" % 
			(output, i), 'w')

		for contig in contigDict:

			temp_num = clusterDict[contig]

			if temp_num == i:

				test_file.write('>'+contig+'\n'+contigDict[contig]+'\n')




# Runs the overall program 
def main(contigFasta, clusterTable, outputPath):

	contig_dict, cluster_dict = create_dictionaires(contigFasta, clusterTable)

	cluster_list = get_cluster_groups(cluster_dict)

	create_output_dir(outputPath)

	create_clustered_fastas(cluster_list, contig_dict, cluster_dict, outputPath)


	

# Initializes at the start of the program
if __name__ == '__main__': 

	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-cf", "--contig_fasta", 
		default="%sall_contigs_1kbto10kb.fasta" % (workdir), 
		type=str, help="Overall combined contig file. \
		The default is set to 'all_contigs_1kbto10kb.fasta'.\n")
	parser.add_argument("-ct", "--cluster_table", 
		default="data/process/clustering_gt1000.csv", 
		type=str, help="Overall cluster table after CONCOCT run. \
		The default is set to 'data/process/clustering_gt1000.csv'.\n")
	parser.add_argument("-o", "--output_dir", 
		default="data/process/clustered_fastas/", 
		type=str, help="Output directory where the clustered fastas will be stored. \
		The default is set to 'data/process/clustered_fastas/'.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.contig_fasta, args.cluster_table, args.output_dir)