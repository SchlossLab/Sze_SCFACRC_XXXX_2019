#!/bin/bash
SEED_DIR=$1

# Keep the first line of File1 and remove the first line of all the others and combine
# This script will combine cv(all) and testing(best) results

echo "seed,mtry,ROC,Sens,Spec,ROCSD,SensSD,SpecSD,model" > $SEED_DIR/mtry_results.csv
grep "Random_Forest" $SEED_DIR/all_hp_results_Random_Forest_*.csv | sed -E 's/.*_Random_Forest_(.*).csv:/\1,/' | sort -k 1,2 -n >> $SEED_DIR/mtry_results.csv

echo "seed,cv_aucs,test_aucs,model" > $SEED_DIR/best_results.csv
grep "Random_Forest" $SEED_DIR/best_hp_results_Random_Forest_*.csv | sed -E 's/.*_Random_Forest_(.*).csv:/\1,/' | sort -k 1,2 -n >> $SEED_DIR/best_results.csv
