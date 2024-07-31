#!/bin/bash

for i in {0..14}
do
	sbatch run_DEBIAS-M.sbatch $i
	echo $i "of 14 done"
done
