#!bash

#Load needed modules
# To find info about the modules use command `module show`
# Needs to be executed with the source command not with bash


# These are from Greg Dick's lab (Niel helped point these out)
# One way to get concoct to work on flux
module use /dept/geology/geomicro/data9/flux/modulefiles
module load geomicro/omics

# Other modules that are needed but can be loaded as normal on flux
module load R/3.3.3
module load sratoolkit/2.8.2-1
module load cutadapt/1.14
module load samtools/1.3.1
module load bowtie2/2.1.0
module load sickle/1.33.6
module load bedtools2/2.20.1
module load megahit/1.0.6
module load bedtools2/2.20.1
module load picard/2.4.1
module load prodigal/2.6.3
	# To run picard tools use java -jar $PICARD_JARS/picard.jar <PROGRAM> <commands>
module load python-anaconda-arc-connect/latest



# Control pathing
export BOWTIE2_INDEXES=/nfs/turbo/schloss-lab/msze/active_projects/metagenome_practice/data/references

