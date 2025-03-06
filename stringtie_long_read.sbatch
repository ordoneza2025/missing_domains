#!/bin/bash

## Example usage
## conda activate Stringtie2
## inFile=. refFile=. outDir=. sbatch stringtie_long_read.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=20:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=stringtie_long_read
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run Stringtie2

echo "Stringtie version: "$(stringtie --version)

pwd; hostname; date

stringtie ${inFile} -p 8 \
-L #long_read mode \
-G ${refFile} \
-A ${outDir}/EF_panel_lr.gene_abund.tab #RSEM quantification of transcripts \
-o ${outDir}/EF_panel_lr.gtf

echo $(date +"[%b %d %H:%M:%S] Done")
