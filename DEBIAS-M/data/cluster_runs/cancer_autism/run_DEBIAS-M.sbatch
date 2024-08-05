#!/bin/bash

#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=after_DEBIAS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1000
#SBATCH --time=20:00

echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
cd $SLURM_SUBMIT_DIR


module load anaconda3/2022.10 # load modules


python DEBIAS-M_semi_leaky.py $1 &> logs/DEBIAS-M_semi_leaky/run_dm/run_DEBIAS-M_out_${1}.log


echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
