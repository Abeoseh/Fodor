#!/bin/bash

if [ -f ./csv_files/AUCs/semi_leaky/semi_leaky_AUCs.csv ]; then
	rm ./csv_files/AUCs/semi_leaky/semi_leaky_AUCs.csv
fi


search_dir=./csv_files/DEBIAS-M_runs/semi_leaky

for entry in "$search_dir/"*debiased*
do
	sbatch run_rf.sbatch $entry 100 $(echo $entry | cut -d"_" -f6- | cut -d"." -f1)

done
echo "done"
