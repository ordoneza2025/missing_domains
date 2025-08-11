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

# Change to output directory for final processing
cd "$outdir"

# Merge all filtered files
cat *.filtered.txt > all_truncations.txt
total_entries=$(wc -l < all_truncations.txt)
echo "Total entries in merged file: $total_entries"

# Collapse duplicate entries
# Key: columns 1-4, Value: columns 5-11
# Combines multiple instances of same key with different values
awk -F'\t' -v OFS='\t' '
{
    key = $1 OFS $2 OFS $3 OFS $4
    value = $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
    if (!seen[key,value]++) {
        data[key] = (key in data ? data[key] OFS : "") value
    }
}
END {
    for (k in data) {
        print k, data[k]
    }
}' all_truncations.txt > collapsed_truncations.txt

collapsed_entries=$(wc -l < collapsed_truncations.txt)

echo "Job completed $(date)"
