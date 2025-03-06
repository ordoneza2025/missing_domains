#!/bin/bash

## Example usage
## short_read_gtf=. long_read_gtf=. hybrid_gtf=. refgtf=. sbatch gffcompare.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=gffcompare
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run GffCompare
pwd; hostname; date

/Users/anor9792/programs/gffcompare/gffcompare ${short_read_gtf} ${long_read_gtf} ${hybrid_gtf} ${refgtf} -M #ignores single exon transcripts

echo $(date +"[%b %d %H:%M:%S] Done")

