#!/bin/bash

for i in {0..2}
do
	sbatch run_DEBIAS-M.sbatch $i
	echo $i "of 2 done"
done
