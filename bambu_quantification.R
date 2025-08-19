#!/usr/bin/env Rscript

################################################################################
# Script: bambu_quantification.R
# Description: Long-read RNA-seq quantification using bambu for human tissue samples

# Dependencies:
#   - R >= 4.0.0
#   - bambu >= 3.0.0
#   - Rsamtools >= 2.8.0

# Input requirements:
#   - BAM files from long-read RNA-seq experiments
#   - Reference genome (FASTA format)
#   - GTF annotation file

# Output:
#   - se.quantOnly_human_FLAIR_counts.txt: Gene/transcript count matrix
#   - se.quantOnly_human_FLAIR_transcript_counts.txt: Transcript-level counts
#   - se.quantOnly_human_FLAIR.gtf: Updated GTF with discovered transcripts
#   - se.quantOnly_human_FLAIRCPM_transcript.txt: CPM values for transcripts (normalized)
#
# Example usage:
# work_dir="/path/to/working/directory" genome_path="/path/to/genome.fa" \
# annotation_path="/path/to/annotation.gtf" bam_dir="/path/to/bam/files" \
# Rscript bambu_quantification.R
#
# Notes:
#   - Filters BAM files for specific organ types: ovary, lung, kidney, pancreas, 
#     pancreatic, colon, mammary, liver
#   - Uses FLAIR-filtered GTF annotation for transcript discovery
#   - Discovery mode disabled for quantification-only analysis
#
################################################################################

## Running bambu quantification for human tissues

# Run job header
cat("Job: bambu_quantification.R\n")
cat("Job started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Load required libraries
suppressPackageStartupMessages({
  library(bambu)
  library(Rsamtools)
})

# Define input parameters and paths
#   work_dir: working directory for analysis
#   genome_path: path to reference genome FASTA file
#   annotation_path: path to GTF annotation file
#   bam_dir: directory containing BAM files
#   output_path: directory for output files
#   output_prefix: prefix for output filenames
work_dir <- Sys.getenv("work_dir", "/scratch/Users/anor9792/bat_decoys/paper/human/missing_domains/domains/filtered/bambu")
genome_path <- Sys.getenv("genome_path", "/Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa")
annotation_path <- Sys.getenv("annotation_path", "/Shares/CL_Shared/db/genomes/hg38/FLAIR_longread_annotation/flair_filter_transcripts.gtf")
bam_dir <- Sys.getenv("bam_dir", "/scratch/Shares/chuong/AI/ENCODE/longread_RNAseq/ENCODE_bams")
output_path <- Sys.getenv("output_path", "quant_only")
output_prefix <- Sys.getenv("output_prefix", "se.quantOnly_human_FLAIR")

# Set working directory
getwd()
setwd(work_dir)

# Load genome and GTF annotation file
cat("Loading genome and annotation files...\n")
genomeSequence <- genome_path
annotation <- prepareAnnotations(annotation_path)

# Define organ pattern for BAM file filtering
#   Pattern matches files containing specific organ names followed by .bam extension
organ_pattern <- "(ovary|lung|kidney|pancreas|pancreatic|colon|mammary|liver).*\\.bam$"

# Get BAM files containing specific organ names
cat("Searching for BAM files with target organ names...\n")
file_paths_panel <- list.files(path = bam_dir, 
                              pattern = organ_pattern, 
                              full.names = TRUE, 
                              ignore.case = TRUE)

cat("Found", length(file_paths_panel), "matching BAM files:\n")
print(file_paths_panel)

# Validate that files exist
if (length(file_paths_panel) == 0) {
  stop("No BAM files found matching the organ pattern. Please check the directory path and file names.")
}

# Read and process BAM files
cat("Creating BamFileList object...\n")
bamFiles_panel <- Rsamtools::BamFileList(file_paths_panel)
head(bamFiles_panel)

# Run bambu quantification
#   discovery: set to FALSE for quantification-only mode
#   reads: BAM files for quantification
#   annotations: prepared GTF annotations
#   genome: reference genome sequence
cat("Starting bambu quantification (discovery mode disabled)...\n")
se.quantOnly <- bambu(reads = bamFiles_panel, 
                     annotations = annotation, 
                     genome = genomeSequence,
                     discovery = FALSE)

# Write output files
#   path: output directory
#   prefix: filename prefix for output files
cat("Writing output files...\n")
writeBambuOutput(se.quantOnly, path = output_path, prefix = output_prefix)

cat("Output files created in directory:", output_path, "\n")
cat("Files prefixed with:", output_prefix, "\n")
cat("Job completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

################################################################################
