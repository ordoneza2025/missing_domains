# Identifying missing domains of alternative bat isoforms 
The purpose of this workflow is to explore alternative isoforms at a transcriptome wide level in two bat species: *Artibeus Jamaicensis* and *Eptesicus Fuscus*. Starting from generated long-read and short-read RNA sequencing, we will assemble a representative transcriptome from a panel of tissues for each of the species. We will then identify instances of alternative isoforms and use protein domain annotation to try to predict something about the function of these isoforms. As most alternative isoforms result in the loss of sequences, we will be focusing on finding which domains would be missing or truncated.    
## Generating a representative transcriptome
## Predicting ORFs from transcripts
## Downloading and re-formatting Uniprot domain files 

These data files are found in the downloadable data for humans in ucsc. 
https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/

1. **UnipDomain.bb** - annotated domains
2. **UnipLocCytopl.bb** - Cytoplasmic domains.
3. **UnipLocExtra.bb** - Extracellular domains.
4. **UnipLocSignal.bb** - signal peptide sequences. 
5. **UnipLocTransMemb.bb** - transmembrane domains.

These files contain both Swissprot (validated) and Trembl (computational annotation) 

The following bash commands will convert .bb files .bed files and re-format them. 
https://github.com/ordoneza2025/missing_domains/blob/main/fromat_Uniprot_domain_files.sh
   
## Identifying missing domains of alternative transcripts - workflow

### 1. Annotating ORFs

The following SLURM script runs BLASTp to annotate the FASTA file of best ORF preidctions for assembled transcripts, using the human swissprot record as the subject database. The output will be BLAST table outfmt 6. 
https://github.com/ordoneza2025/missing_domains/blob/main/blastp.sbatch

### 2. Finding truncations in bat transcripts

The following script parses through the human swissprot and the bat best ORFs FASTA files. It uses the BLASTp output to match sequences, which are then aligned by MUSCLE. Any alignment gaps larger than 5 amino acids are flagged as truncations. Coordinates for the gap regions correspondant to human amino acid sequences are written into a BED format file. 
https://github.com/ordoneza2025/missing_domains/blob/main/identifying_truncations.sbatch

### 3. Assigning truncations to protein domains

The following script uses bedtools to intersect the truncations identifies in STEP2 with protein domains, based on amino acid coordinates. The commands will be ran for each of the 5 domain files described above.   








