#!/bin/bash

if [ -f ./csv_files/AUCs.csv ]; then
	rm ./csv_files/AUCs.csv
fi

for i in {1..15}
do
	Rscript random_forest.R $i 15 &> logs/run_rf_out_$i.log
	echo $i "of 3 done"
done