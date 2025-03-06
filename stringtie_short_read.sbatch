#!/bin/bash

## Example usage
## conda activate Stringtie2
## inFile=. refFile=. outDir=. sbatch stringtie_short_read.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=20:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=stringtie_short_read
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run Stringtie2

echo "Stringtie version: "$(stringtie --version)
pwd; hostname; date

stringtie ${inFile} -p 8 #threads \
-A ${outDir}/AJ_all_sr.gene_abund.tab #RSEM transcript quantification \
--rf #library strandedness \
-G ${refFile} \
-o ${outDir}/AJ_panel_sr.gtf

echo $(date +"[%b %d %H:%M:%S] Done")
