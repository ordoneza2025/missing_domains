#!/bin/bash

#SBATCH --job-name=missing_domains
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=/scratch/Users/anor9792/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/anor9792/eofiles/%x.%j.err

#Activate Biopython environment
source /Users/anor9792/biopython/biopython/bin/activate

# Execute the Python script
python missing_domains_muscle3.py


echo $(date +"[%b %d %H:%M:%S] Done")



