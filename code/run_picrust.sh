#!bash

#Set local variables
PICRUST=~/picrust/scripts
WORKDIR=data/process

# create a useable biom formated table for picrust
# help from https://github.com/rprops/PICRUSt_from_mothur

# Create an empty nsti file for later use
touch $WORKDIR/nsti.txt

# Comments on whether this is an acceptable biom format
biom validate-table -i $WORKDIR/gg_final.0.03.biom
# Convert to an acceptable biom format for picrust
biom convert --table-type="OTU table" -i $WORKDIR/gg_final.0.03.biom -o $WORKDIR/temp_OTU_table.txt --to-tsv --header-key taxonomy
biom convert -i $WORKDIR/temp_OTU_table.txt -o $WORKDIR/gg_corr_OTU_table.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy
# Double check that it is now the right format
biom validate-table -i $WORKDIR/gg_corr_OTU_table.biom
# Run the picrust program to make predictions 
$PICRUST/normalize_by_copy_number.py -i $WORKDIR/gg_corr_OTU_table.biom -o $WORKDIR/gg_corr_normalized_otus.biom
$PICRUST/predict_metagenomes.py -i $WORKDIR/gg_corr_normalized_otus.biom -o $WORKDIR/predicted_metagenomes.biom -a $WORKDIR/nsti.txt
		# using -f tag creates a tab delimited table

#biom convert --table-type="OTU table" -i normalized_otus.biom -o normalized_otus_table.txt --header-key taxonomy --to-tsv
#biom convert -i normalized_otus.biom -o normalized_otus_table.txt  --table-type "otu table"
