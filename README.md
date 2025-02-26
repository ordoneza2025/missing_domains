# Identifying missing domains of alternative bat isoforms 
The purpose of this workflow is to explore alternative isoforms at a transcriptome wide level in two bat species: *Artibeus Jamaicensis* and *Eptesicus Fuscus*. Starting from generated long-read and short-read RNA sequencing, we will assemble a representative transcriptome from a panel of tissues for each of the species. We will then identify instances of alternative isoforms and use protein domain annotation to try to predict something about the function of these isoforms. As most alternative isoforms result in the loss of sequences, we will be focusing on finding which domains would be missing or truncated.    
## Generating a representative transcriptome
## Predicting ORFs from transcripts
## Identifying missing domains of alternative transcripts - workflow

### 1. Annotating ORFs

The following SLURM script runs BLASTp to annotate the FASTA file of best ORF preidctions for assembled transcripts, using the human swissprot record as the subject database. The output will be BLAST table outfmt 6. 
https://github.com/ordoneza2025/missing_domains/blob/main/blastp.sbatch

### 2. Finding truncations in bat transcripts

The following PYTHON script parses through the human swissprot and the bat best ORFs FASTA files. It uses the BLASTp output to match sequences, which are then aligned by MUSCLE. Any alignment gaps larger than 5 amino acids are flagged as truncations. Coordinates for the gap regions correspondant to human amino acid sequences are written into a BED format file. The script is submitted via SLURM.






