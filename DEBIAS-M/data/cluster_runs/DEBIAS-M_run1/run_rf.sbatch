#!/bin/bash

#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=after_DEBIAS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000
#SBATCH --time=0:20:00

echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
cd $SLURM_SUBMIT_DIR


module load R # load modules


Rscript random_forest.R $1 $2 &> logs/run_rf_out_${1}.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
