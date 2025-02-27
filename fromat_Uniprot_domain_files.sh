#!/bin/bash

#Example commands to convert and re-format .bb uniprot domain files to .bed

#bigBedtoBed command in part of the UCSC Genome Browser Utilities

#STEP1: Loop through all .bb files in the current directory
for file in *.bb; do
    # Check if any .bb files exist
    [ -e "$file" ] || continue
    
    # Define output .bed file name
    output="${file%.bb}.bed"
    
    # Convert bigBed to BED
    bigBedToBed "$file" "$output"
    
    echo "Converted $file to $output"
done

## STEP2: choose relevant columns 
awk -F'\t' '{print $28 "\t" $27 "\t" $23}' UnipDomain.bed > domains.txt

## STEP3: simplify column 2 to just amino acid ranges
awk -F'\t' '{
    if (match($2, /amino (acid|acids) ([0-9]+)-([0-9]+)/, arr)) {
        print $1 "\t" arr[2] "-" arr[3] "\t" $3
    } else if (match($2, /amino (acid|acids) ([0-9]+)/, arr)) {
        print $1 "\t" arr[2] "\t" $3
    }
}' domains.txt > domains.2.txt

## STEP #4: convert amino acid ranges to BED format. 
sed -E 's/^([^ \t]+)[ \t]+([0-9]+)-([0-9]+)/\1\t\2\t\3/' domains.2.txt > domains.3.bed
sort domains.3.bed > domains.4.bed

# STEP5: AT this point I noticed that some of the domains are just a single aa, so there is no range, meaning incorrect number of columns for BED format. Keep only the entries with a valid range. The invalid entries don't add important infromation, most are just described as lumenal and there are only 45 of them.  
awk -F'\t' 'NF==4' domains.4.bed  > valid_domains.4.bed

##Followed STEPS 2-4 for UnipLocSignal.bed file , all lines are a range of aa so skipped step 5.

##For UnipLocCytopl.bed, follow STEP2, with output cyto.txt. STEP3 and STEP4 will be slightly different.

##STEP 3for cyto.txt
awk -F'\t' '{
    if (match($3, /amino (acid|acids) ([0-9]+)-([0-9]+)/, arr)) {
        print $1 "\t" arr[2] "-" arr[3] "\t" $2
    } else if (match($3, /amino (acid|acids) ([0-9]+)/, arr)) {
        print $1 "\t" arr[2] "\t" $2
    }
}' cyto.txt > cyto.2.txt

##STEP4 for cyto.2.txt
sed -E 's/^([^ \t]+)[ \t]+([0-9]+)-([0-9]+)/\1\t\2\t\3/' cyto.2.txt | sort > cyto.3.bed

##UnipLocTransMem.bed an UnipLocExtra.bed were processed with the sames steps as UnipLocCytopl.bed

##The last step is to add the protein names to each bed file. Starting from the dowloaded swissprot.fasta record
awk -F '[| ]+' '/^>/ {split($3, name, "_"); print $2 "\t" name[1]}' swissprot.fasta > protein_names.txt

# Define the mapping file
map_file="protein_names.txt"

# Loop through all .bed files in the directory
for bedfile in *.bed; do
    # Skip the output files to avoid reprocessing
    [[ "$bedfile" == *.protnames.bed ]] && continue

    # Define output file name
    output="${bedfile%.bed}.protnames.bed"

    # Run awk command
    awk 'BEGIN { OFS="\t" } NR==FNR { map[$1]=$2; next } { print $0, map[$1] }' "$map_file" "$bedfile" > "$output"

    echo "Processed $bedfile -> $output"
done
