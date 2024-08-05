#!/bin/bash

if [ -f ./csv_files/AUCs/semi_leaky_autism/autism_semi_leaky_AUCs.csv ]; then
	rm ./csv_files/AUCs/semi_leaky_autism/autism_semi_leaky_AUCs.csv
fi


search_dir=./csv_files/DEBIAS-M_runs/semi_leaky_autism

for entry in "$search_dir/"*debiased*
do
	sbatch run_rf_post_DEBIAS.sbatch $entry 100 $(echo $entry | cut -d"_" -f8- | cut -d"." -f1)
	#(echo $entry | cut -d"_" -f8- | cut -d"." -f1)
	#echo $entry
done
echo "done"
