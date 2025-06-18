#!/bin/bash

## Example usage
## input_bed=. bed_files_dir=. output_dir=. sbatch uniprot_domains_intersect.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=uniprot_domains_intersect
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run job header
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on $(hostname)
echo Job started $(date)
echo Using $SLURM_NTASKS processors and $SLURM_NNODES nodes

#Load necessary modules
module load bedtools

# Loop through all .bed files in the directory
for bed_file in "${bed_files_dir}"/*.bed; do
    # Extract the base filename without the path and remove the .bed extension
    base_name=$(basename "$bed_file" .bed)
    
    # Define the output file with .txt extension
    output_file="${output_dir}/truncations_${base_name}.txt"
    
    # Perform the intersection and append to output
    bedtools intersect -a "${input_bed}" -b "$bed_file" -wo >> "$output_file"
    
    # Process the file to add a new first column with the TCONS_# ID (without .p#)
    awk -F'\t' '{
        match($4, /(TCONS_[0-9]+)/, id);
        print id[1], $0;
    }' OFS='\t' "$output_file" > "${output_file}.tmp" && mv "${output_file}.tmp" "$output_file"
    
    # Notify completion for this file
    echo "Intersection and processing completed for ${base_name}"
done

# Final message
echo "All intersections are complete."


#| awk '{print $4, $1, $2, $3, $8, $6, $7}' OFS='\t'