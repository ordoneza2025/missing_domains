#!/bin/bash

################################################################################
# Script: CPM_filtering.sbatch
# Description: Filtering transcripts based on CPM from quantified tissue panel

# Dependencies:


# Input requirements:
#   - CPM.txt file output from bambu


# Output:
#   - selected_transcript_rows_with_header.txt: All transcripts with header
#   - transcript_max_cpm_by_tissue.txt: max CPM from correspondant tissue for each transcript
#   - transcript_CP5.txt: filtered transcripts with a threshold of 5 CPM

# Example usage:
# input_file="data.txt" cpm_threshold=5 sbatch transcript_cpm_analysis.sbatch

# Notes:
#   - Removes .bam suffix from tissue column headers
#   - Identifies tissue with maximum CPM for each transcript
#   - Filters results by minimum CPM threshold

################################################################################

## Example usage
## input_file="se.quantOnly_human_FLAIRCPM_transcript.txt" cpm_threshold=5 sbatch transcript_cpm_analysis.sbatch

## General Settings
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anor9792@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --partition=short

## Job name and output
#SBATCH --job-name=transcript_cpm_analysis
#SBATCH --output=/scratch/Users/anor9792/eofiles/%x.%j.out
#SBATCH --error=/scratch/Users/anor9792/eofiles/%x.%j.err

# Run job header
echo "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
echo "Running on $(hostname)"
echo "Job started $(date)"
echo "Using $SLURM_NTASKS processors and $SLURM_NNODES nodes"

# Define input parameters
#   input_file: CPM transcript file with tissue expression data
#   cpm_threshold: minimum CPM value for filtering (default: 5)
#   temp_file1: intermediate file with header and all transcripts
#   temp_file2: file with max CPM per transcript
#   output_file: final filtered file
input_file="${input_file:-se.quantOnly_human_FLAIRCPM_transcript.txt}"
cpm_threshold="${cpm_threshold:-5}"
temp_file1="selected_transcript_rows_with_header.txt"
temp_file2="transcript_max_cpm_by_tissue.txt"
output_file="transcript_CPM${cpm_threshold}.txt"

# Validate input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file '$input_file' not found"
    exit 1
fi

echo "Processing input file: $input_file"
echo "CPM threshold: $cpm_threshold"
echo ""

# Step 1: Extract all transcript rows with header
#   Copies all rows from input file (header + all transcripts)
echo "Step 1: Extracting transcript rows with header..."
awk 'NR == 1 { print $0; next }  # Print the header first
     { print $0 }' "$input_file" > "$temp_file1"

transcript_count=$(tail -n +2 "$temp_file1" | wc -l)
echo "Extracted $transcript_count transcripts"

# Step 2: Find maximum CPM and corresponding tissue for each transcript
#   Analyzes columns 3-NF to find tissue with highest expression
#   Removes .bam suffix from tissue names for cleaner output
echo "Step 2: Finding maximum CPM and tissue for each transcript..."
awk -F'\t' -v OFS='\t' '
NR==1 {
  # build an array of tissue names from columns 3…NF
  for(i=3;i<=NF;i++){
    h[i]=$i
    sub(/\.bam$/,"",h[i])
  }
  next
}
{
  # keep the first two columns (transcript info)
  transcript_info = $1 OFS $2
  # initialize with column 3
  max = $3
  tissue = h[3]
  # scan columns 4 to NF to find maximum
  for(i=4;i<=NF;i++){
    if($i+0 > max+0){
      max = $i
      tissue = h[i]
    }
  }
  print transcript_info, max, tissue
}' "$temp_file1" > "$temp_file2"

max_analysis_count=$(wc -l < "$temp_file2")
echo "Analyzed maximum CPM for $max_analysis_count transcripts"

# Step 3: Filter transcripts by CPM threshold
#   Keeps only transcripts with maximum CPM >= threshold
echo "Step 3: Filtering transcripts with CPM >= $cpm_threshold..."
awk -F'\t' -v OFS='\t' -v threshold="$cpm_threshold" '$3+0 >= threshold' "$temp_file2" > "$output_file"

filtered_count=$(wc -l < "$output_file")
echo "Filtered to $filtered_count transcripts meeting CPM threshold"

echo "Job completed $(date)"
