#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200G
#SBATCH --partition=Orion

cd $SLURM_SUBMIT_DIR

module load kraken2

# pass the file name as an arg
FILE1=$1
# get the reverse file name
FILE2=${FILE1/_R1_/_R2_}

echo "Preparing to process paired files: $FILE1 and $FILE2."
du -sh $FILE1
du -sh $FILE2

# set up output path
OUTDIR=/projects/afodor_research2/ieclabau/bolivia_2024/kraken2/output
BASENAME=$(basename $FILE1)
OUTNAME=${BASENAME/_R1_001.fastq.gz/.kraken}
OUT=$OUTDIR/$OUTNAME

REPORTNAME=${BASENAME/_R1_001.fastq.gz/.kreport}
REPORT=$OUTDIR/$REPORTNAME
echo "Output will be saved to: $OUT"
echo "Report will be saved to: $REPORT"


DBNAME=/projects/afodor_research2/databases/kraken2_Standard_12-2024
echo "Database is: $DBNAME"

kraken2 \
   --threads 8 \
   --gzip-compressed \
   --db $DBNAME \
   --report $REPORT \
   --paired \
   $FILE1 $FILE2 > $OUT

echo "Done."
