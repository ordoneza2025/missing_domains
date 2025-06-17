
## Example usage
## blastp_output=. best_orfs_fasta=. swissprot_fasta=.  Submit job through slurm by running sbatch identifying_truncations.sbatch

#Import programs requires to handle sequence info and carry out alignment
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import subprocess
import os

# Parse the BLASTp output file to extract Bat protein and SwissProt IDs
with open(blastp_output, "r") as blastp_file, open(parsed_blastp_file, "w") as outfile:
    for line in blastp_file:
         fields = line.strip().split("\t")
         bat_protein = fields[0]
         #swissprot_id = fields[1]
         swissprot_id = fields[1].split("|")[1]
         outfile.write(f"{bat_protein}\t{swissprot_id}\n")


# Step 2: Parse human and bat sequences

best_orfs_fasta = "AJ_bestorfs_pep.fasta"  # Bat protein sequences
swissprot_fasta = "swissprot.fasta"  # Human protein sequences
parsed_blastp_file = "parsed_blastp.txt"  # Input from step 1
sequences_file = "sequences.txt"  # Output file for sequences

# Load Bat and Human sequences from FASTA files
bat_sequences = {record.id: record.seq for record in SeqIO.parse(best_orfs_fasta, "fasta")}
human_sequences = {record.id.split('|')[1]: record.seq for record in SeqIO.parse(swissprot_fasta, "fasta")}

# Get sequences for each Bat and SwissProt ID pair and write to file
with open(parsed_blastp_file, "r") as infile, open(sequences_file, "w") as outfile:
    for line in infile:
        bat_protein, swissprot_id = line.strip().split("\t")
        if bat_protein in bat_sequences and swissprot_id in human_sequences:
            bat_sequence = bat_sequences[bat_protein]
            human_sequence = human_sequences[swissprot_id]
            outfile.write(f"{bat_protein}\t{swissprot_id}\t{bat_sequence}\t{human_sequence}\n")

# Step 3: Running muscle3 alignment, default output 

# Define file paths
alignment_file = "muscle3.alignment.txt"
sequences_file = "sequences.txt"
checkpoint_file = "checkpoint.txt"

# Load the checkpoint
try:
    with open(checkpoint_file, "r") as ckpt_file:
        last_processed_line = int(ckpt_file.read().strip())
except FileNotFoundError:
    last_processed_line = 0  # Start from the beginning if no checkpoint exists

print(f"Starting from line {last_processed_line + 1}...")

# Alignment process
with open(sequences_file, "r") as infile, open(alignment_file, "a") as outfile:
    for line_number, line in enumerate(infile, 1):  # Enumerate for line tracking
        if line_number <= last_processed_line:
            continue  # Skip already processed lines

        print(f"Processing line {line_number}: {line.strip()}")
        try:
            bat_protein, swissprot_id, bat_sequence, human_sequence = line.strip().split("\t")

            # Write sequences to a temporary FASTA file
            temp_fasta = f"temp_alignment_{line_number}.fasta"
            with open(temp_fasta, "w") as temp_file:
                temp_file.write(f">Bat_{bat_protein}\n{bat_sequence.replace('*', '')}\n")
                temp_file.write(f">Human_{swissprot_id}\n{human_sequence}\n")

            # Define MUSCLE output file for the alignment
            temp_alignment_output = f"temp_alignment_{line_number}.txt"

            # Run MUSCLE for alignment
            try:
                result = subprocess.run(
                    ["muscle", "-in", temp_fasta, "-out", temp_alignment_output],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                )
                print(f"MUSCLE output:\n{result.stdout}")
            except subprocess.CalledProcessError as e:
                print(f"Error running MUSCLE for line {line_number}: {e.stderr}")
                continue  # Skip to the next sequence pair if MUSCLE fails

            # Read the MUSCLE alignment result
            with open(temp_alignment_output, "r") as alignment_file:
                alignment_data = alignment_file.read()

            # Parse alignment result and write to the output file
            outfile.write(f"{bat_protein}\t{swissprot_id}\n{alignment_data}\n")
            print(f"Line {line_number} processed successfully.")

            # Clean up temporary files
            os.remove(temp_fasta)
            os.remove(temp_alignment_output)

        except Exception as e:
            print(f"Error processing line {line_number}: {e}")
            continue

        # Update the checkpoint after each successfully processed line
        with open(checkpoint_file, "w") as ckpt_file:
            ckpt_file.write(str(line_number))

#Step 4: re-format output so that the protein ids and the aligned sequences are seperated by tab deliminator

input_file = "muscle3.alignment.txt"
output_file = "rearranged_musc3_output.txt"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    current_pair = {}
    current_sequences = {"Bat": "", "Human": ""}
    
    for line in infile:
        line = line.strip()
        if line.startswith(">Bat_"):
            current_type = "Bat"
        elif line.startswith(">Human_"):
            current_type = "Human"
        elif line.startswith("TCONS_"):
            # Write the previous pair if it exists
            if current_pair:
                aligned = f"{current_pair['id']}\t{current_pair['ref_id']}\t{current_sequences['Bat']}\t{current_sequences['Human']}"
                outfile.write(aligned + "\n")
                current_sequences = {"Bat": "", "Human": ""}
            # Start a new pair
            parts = line.split()
            current_pair = {"id": parts[0], "ref_id": parts[1]}
        else:
            # Append the sequence
            current_sequences[current_type] += line
    
    # Write the last pair
    if current_pair:
        aligned = f"{current_pair['id']}\t{current_pair['ref_id']}\t{current_sequences['Bat']}\t{current_sequences['Human']}"
        outfile.write(aligned + "\n")


#Step 5: create bed file with truncated domain ranges (gaps larger than 10 amino acids)

alignment_file = "rearranged_musc3_output.txt"  # Input from step 3
output_file = "missing_regions_human.bed"  # Final output in BED format


# Open the alignment file for reading and output file for writing
with open(alignment_file, "r") as infile, open(output_file, "w") as outfile:
    # Iterate through each line in the input file
    for line in infile:
        # Parse the input line into its respective fields
        bat_protein, swissprot_id, aligned_bat, aligned_human = line.strip().split("\t")
        
        gap_start = None  # Initialize gap start as None
        for i in range(len(aligned_bat)):
            if aligned_bat[i] == "-" and gap_start is None:
                # Mark the start of a new gap
                gap_start = i
            elif (aligned_bat[i] != "-" or i == len(aligned_bat) - 1) and gap_start is not None:
                # Handle the end of the gap or the end of the alignment
                gap_end = i if aligned_bat[i] != "-" else i + 1
                gap_length = gap_end - gap_start
                
                # Process gaps larger than 10
                if gap_length > 10:
                    # Extract the corresponding region from the aligned human sequence
                    human_gap_region = aligned_human[gap_start:gap_end]
                    
                    # Calculate the indices in the human sequence (excluding '-')
                    human_gap_start = len(aligned_human[:gap_start].replace("-", ""))
                    human_gap_end = human_gap_start + len(human_gap_region.replace("-", ""))
                    
                    # Write the simplified output to the file
                    outfile.write(
                        f"{swissprot_id}\t{human_gap_start}\t{human_gap_end}\t{bat_protein}\n"
                    )
                
                # Reset gap start for the next potential gap
                gap_start = None


#Step 6: load program to merge adjacent gaps into one.

from collections import defaultdict

input_bed = "missing_regions_human.bed"
merged_bed = "missing_regions_human_merged.bed"

#Step 7: Merge nearby gaps (≤15 aa apart) per bat isoform
gaps_by_protein = defaultdict(list)

# Read original bed entries
with open(input_bed, "r") as infile:
    for line in infile:
        ref_id, start, end, bat_protein = line.strip().split("\t")
        gaps_by_protein[bat_protein].append((ref_id, int(start), int(end)))

# Merge gaps ≤15 aa apart
with open(merged_bed, "w") as outfile:
    for bat_protein, gaps in gaps_by_protein.items():
        # Sort gaps by start coordinate
        sorted_gaps = sorted(gaps, key=lambda x: x[1])
        merged = []
        
        current_ref, current_start, current_end = sorted_gaps[0]
        
        for ref_id, start, end in sorted_gaps[1:]:
            if start - current_end <= 15:
                # Extend the current gap
                current_end = max(current_end, end)
            else:
                # Finalize the previous gap and start a new one
                merged.append((current_ref, current_start, current_end, bat_protein))
                current_ref, current_start, current_end = ref_id, start, end
        
        # Add the last gap
        merged.append((current_ref, current_start, current_end, bat_protein))
        
        # Write to file
        for ref_id, start, end, bp in merged:
            outfile.write(f"{ref_id}\t{start}\t{end}\t{bp}\n")
