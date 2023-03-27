#!/bin/bash

#this script downloads samples from sra -> trims adapters and does quality control -> aligns and quantifies paired end rna seq samples
filename='samples_PE.txt' #txt file with list of all sample names (ie SRR number) with \n between them
echo Start
while read p; do #loop through each sample
    prefetch $p; 
    fastq-dump --split-files *.sra; #download it from the SRA
    trim_galore --fastqc --gzip --paired --cores 8 --output_dir . "${p}_1.fastq" "${p}_2.fastq"; #trim adapters and low quality nucleotides
    mv *.html ./fastqc_files/PE; #move html file containing quality info about sample into different folder 
    salmon quant -i ./translatome_k_25_idx -p 20 --validateMappings --writeUnmappedNames -l A -1 "${p}_1_val_1.fq.gz" -2 "${p}_2_val_2.fq.gz" -o "./counts/PE/${p}"; #align and quantify using Salmon
done < $filename
