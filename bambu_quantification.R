#!/usr/bin/env Rscript

################################################################################
# Script: bambu_quantification.R
# Description: Long-read RNA-seq quantification using bambu for human tissue samples
# Author: [Your Name]
# Date: [Date]
# Version: 1.0.0
#
# Dependencies:
#   - R >= 4.0.0
#   - bambu >= 3.0.0
#   - Rsamtools >= 2.8.0
#
# Input requirements:
#   - BAM files from long-read RNA-seq experiments
#   - Reference genome (FASTA format)
#   - GTF annotation file
#
# Output:
#   - LR_human_panel_counts.txt: Gene/transcript count matrix
#   - LR_human_panel_transcript_counts.txt: Transcript-level counts
#   - LR.Human.panel.gtf: Updated GTF with discovered transcripts
#
# Example usage:
# Rscript bambu_quantification.R
#
# Notes:
#   - Filters BAM files for specific organ types: ovary, lung, kidney, pancreas, 
#     pancreatic, colon, mammary, liver
#   - Uses FLAIR-filtered GTF annotation for transcript discovery
#   - Automatically discovers novel transcripts during quantification
#
################################################################################

## Running bambu quantification for human tissues

# Load required libraries
suppressPackageStartupMessages({
  library(bambu)
  library(Rsamtools)
})

# Set working directory
getwd()
setwd("/scratch/Users/anor9792/bat_decoys/paper/human/missing_domains/domains/filtered/bambu")

# Load genome and GTF annotation file
genomeSequence <- "/Shares/CL_Shared/db/genomes/hg38/fa/hg38.main.fa"
annotation <- prepareAnnotations("/Shares/CL_Shared/db/genomes/hg38/FLAIR_longread_annotation/flair_filter_transcripts.gtf")

# Specify the directory containing the BAM files
bam_directory <- "/scratch/Shares/chuong/AI/ENCODE/longread_RNAseq/ENCODE_bams"

# Get BAM files containing specific organ names
organ_pattern <- "(ovary|lung|kidney|pancreas|pancreatic|colon|mammary|liver).*\\.bam$"
file_paths_panel <- list.files(path = bam_directory, 
                              pattern = organ_pattern, 
                              full.names = TRUE, 
                              ignore.case = TRUE)

print(file_paths_panel)

# Validate that files exist
if (length(file_paths_panel) == 0) {
  stop("No BAM files found matching the organ pattern. Please check the directory path and file names.")
}

# Read and process BAM files
bamFiles_panel <- Rsamtools::BamFileList(file_paths_panel)
head(bamFiles_panel)

# Run bambu quantification
LR.Human.panel <- bambu(reads = bamFiles_panel, 
                       annotations = annotation, 
                       genome = genomeSequence)

# Write output files
writeBambuOutput(LR.Human.panel, prefix = "LR_human_panel")
writeToGTF(rowRanges(LR.Human.panel), file = "LR.Human.panel.gtf")

################################################################################
