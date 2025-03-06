#!/bin/bash

## Example usage
## conda activate Stringtie2
## inFile=. lrFile=. srFile=. outDir=. outfilename=. sbatch stringtie_hybrid.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=20:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=stringtie_hybrid
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run Stringtie2

echo "Stringtie version: "$(stringtie --version)
pwd; hostname; date

stringtie -p 8 #threads\
--mix #hybrid_mode \
-o ${outDir}/${outfilename} \
-G ${refFile} \
${srFile} \
${lrFile} \
-A ${outDir}/gene_abund_spleen.tab #RSEM quantification of transcripts
