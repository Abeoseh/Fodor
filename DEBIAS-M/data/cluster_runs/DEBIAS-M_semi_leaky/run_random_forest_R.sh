#!/bin/bash

if [ -f ./csv_files/AUCs.csv ]; then
	rm ./csv_files/AUCs.csv
fi

for i in {1..15}
do
	sbatch run_rf.sbatch $i 100
	echo $i "of 15 done"
done
