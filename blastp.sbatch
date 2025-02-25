#!/bin/bash

## Example usage
## query=. subject_DB=. outDir=. sbatch blastp.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32gb
#SBATCH --time=15:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=balstp
#SBATCH --output=/scratch/Users/anor9792/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/anor9792/eofiles/%x.%j.err

# Set variables
evalue=1e-50
outName=AJ_bestorfs_swissprot.txt
numThreads=16

# Load Module
module load ncbi-blast/2.7.1

# Run blastp
pwd; hostname; date

blastp -query ${query} -db ${subject_DB} -out ${outDir}/${outName}.txt \
-evalue ${evalue} \
-num_threads ${numThreads} -max_target_seqs 1 -outfmt 6

echo $(date +"[%b %d %H:%M:%S] Done")
