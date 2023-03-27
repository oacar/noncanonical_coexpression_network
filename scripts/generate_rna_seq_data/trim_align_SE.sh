#!/bin/bash
#this script downloads samples from sra -> trims adapters and does quality control -> aligns and quantifies single end rna seq samples
filename='samples_SE.txt'  #txt file with list of all sample names (ie SRR number) with \n between them
echo Start
while read p; do #loop through each sample
    fasterq-dump $p; #download it from the SRA
    trim_galore --fastqc --gzip --cores 8 --output_dir . "${p}.fastq"; #trim adapters and low quality nucleotides
    mv *.html ./fastqc_files/SE; #move html file containing quality info about sample into different folder 
    salmon quant -i /home/aar75/rna_seq/Salmon_20221011/k_25_idx -p 20 --validateMappings --writeUnmappedNames -l A -r "${p}_trimmed.fq.gz" -o "./counts/SE/${p}"; #align and quantify using Salmon

done < $filename
~
~
~
~
~
~
~
~
~