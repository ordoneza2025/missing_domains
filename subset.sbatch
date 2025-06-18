#!/bin/bash

## Example usage
## indir=. outdir=. sbatch subset_domains_intersect.sbatch

##General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --partition=short

##Job name and output
#SBATCH --job-name=subset_domains_intersect
#SBATCH --output=/scratch/Users/%u/data/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/%u/data/eofiles/%x.%j.err

# Run job header
echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on $(hostname)
echo Job started $(date)
echo Using $SLURM_NTASKS processors and $SLURM_NNODES nodes

# Define input and output directories
#   indir : directory containing domain intersection files
#   outdir: directory where filtered files will be written
# Make sure outdir exists
mkdir -p "$outdir"

# Declare an associative array for thresholds
#   keys correspond to domain types; values are truncation thresholds
declare -A thresholds=(
    [cyto]=50
    [extra]=50
    [signal]=20
    [tm]=15
    [valid_domains]=50
)

# Loop through each domain type
for domain in "${!thresholds[@]}"; do
    # Identify matching file
    file=$(ls "$indir"/truncations_"$domain".*.protnames.txt 2>/dev/null)
    
    # Skip if file doesn't exist
    if [[ -z "$file" ]]; then
        echo "No file found for domain: $domain"
        continue
    fi

    # Extract base filename
    basefile=$(basename "$file")

    # Set output file path
    outfile="$outdir/${basefile%.txt}.filtered.txt"

    # Apply threshold filter
    awk -v threshold="${thresholds[$domain]}" '$NF >= threshold' "$file" > "$outfile"

    # Notify completion for this file
    echo "Filtered $file with threshold ${thresholds[$domain]} -> $outfile"
done

# Final message
echo "All domain files have been filtered."
