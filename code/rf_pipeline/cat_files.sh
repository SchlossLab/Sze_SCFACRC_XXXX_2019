#!/bin/bash
SEARCH_DIR=data/temp
FINAL_DIR=data/process
# Keep the first line of File1 and remove the first line of all the others and combine
# This script will combine cv(all) and testing(best) results

for model in "Random_Forest"
do
  	head -1 $SEARCH_DIR/all_hp_results_"$model"_1.csv  > $SEARCH_DIR/combined_all_hp_results_"$model".csv; tail -n +2 -$

    head -1 $SEARCH_DIR/best_hp_results_"$model"_1.csv  > $SEARCH_DIR/combined_best_hp_results_"$model".csv; tail -n +2$

    mv $SEARCH_DIR/combined_all_hp_results_"$model".csv $FINAL_DIR/combined_all_hp_results_"$model".csv
    mv $SEARCH_DIR/combined_best_hp_results_"$model".csv $FINAL_DIR/combined_best_hp_results_"$model".csv
done

rm $SEARCH_DIR/all_hp_*
rm $SEARCH_DIR/best_hp_*
