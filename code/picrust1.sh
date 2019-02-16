#!bash

# module load python-anaconda2/latest

# module load R/3.5.1

#Set local variables
REF=data/references
WORKDIR=data/picrust1
mkdir -p $WORKDIR

# create a useable biom formated table for picrust
# help from https://github.com/rprops/PICRUSt_from_mothur

BIOM=$WORKDIR/crc.0.03.biom

R -e "source('code/picrust1_utilities.R');clean_tax()"

cp data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared $WORKDIR/crc.shared

mothur "#list.otus(shared=$WORKDIR/crc.shared);get.otus(constaxonomy=$WORKDIR/crc.cons.taxonomy, accnos=$WORKDIR/crc.0.03.otulabels);make.biom(shared=current, constaxonomy=current, label=0.03, reftaxonomy=$REF/gg_13_5_99.gg.tax, picrust=$REF/97_otu_map.txt)"


source activate picrust1


# Comments on whether this is an acceptable biom format
biom validate-table -i $BIOM

# Convert to an acceptable biom format for picrust
biom convert --table-type="OTU table" -i $BIOM -o $WORKDIR/temp_OTU_table.txt --to-tsv --header-key taxonomy
biom convert -i $WORKDIR/temp_OTU_table.txt -o $BIOM --table-type="OTU table" --to-json --process-obs-metadata taxonomy

# Double check that it is now the right format
biom validate-table -i $BIOM

# Run the picrust v1 program to make predictions
normalize_by_copy_number.py -i $BIOM -o $WORKDIR/crc.normalized.biom
predict_metagenomes.py -f -i $WORKDIR/crc.normalized.biom -o $WORKDIR/crc.metagenomes.tsv -a $WORKDIR/crc.nsti

source deactivate

R -e "source('code/picrust1_utilities.R');get_shared_annotation()"
