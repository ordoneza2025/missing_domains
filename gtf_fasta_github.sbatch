#!/bin/bash

## Example usage
## gtfFile=. fastaGenome=. output=. sbatch gtf_fasta.sbatch

##General Settings
#SBATCH --job-name=GTF_to_FASTA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=12gb
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=/scratch/Users/anor9792/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/anor9792/eofiles/%x.%j.err


# Load Module
module load cufflinks 
module load samtools 


# Index fasta genome file
samtools faidx ${fasta}

#Run gffread from cufflinks 
gffread -w ${output} -g ${fastaGenome} ${gtfFile}
