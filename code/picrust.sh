#!bash

#Set local variables
REF=data/references
WORKDIR=data/picrust
mkdir -p $WORKDIR

# create a useable biom formated table for picrust
# help from https://github.com/rprops/PICRUSt_from_mothur

# Create an empty nsti file for later use
#touch $WORKDIR/nsti.txt

BIOM=data/picrust/crc.0.03.biom

# make biom format from subsampled shared file
Rscript code/picrust_align_metadata.R #this uses the subsampled shared file
mothur "#make.biom(shared=data/picrust/crc.shared, label=0.03, reftaxonomy=$REF/gg_13_5_99.gg.tax, constaxonomy=data/picrust/crc.cons.taxonomy, metadata=data/picrust/crc.metadata, picrust=$REF/97_otu_map.txt)"

# Comments on whether this is an acceptable biom format
biom validate-table -i $BIOM

# Convert to an acceptable biom format for picrust
biom convert --table-type="OTU table" -i $BIOM -o $WORKDIR/temp_OTU_table.txt --to-tsv --header-key taxonomy
biom convert -i $WORKDIR/temp_OTU_table.txt -o $BIOM --table-type="OTU table" --to-json --process-obs-metadata taxonomy

# Double check that it is now the right format
biom validate-table -i $BIOM

# Run the picrust program to make predictions
normalize_by_copy_number.py -i $BIOM -o $WORKDIR/crc.normalized.biom
predict_metagenomes.py -f -i $WORKDIR/crc.normalized.biom -o $WORKDIR/crc.metagenomes.tsv -a $WORKDIR/crc.nsti
