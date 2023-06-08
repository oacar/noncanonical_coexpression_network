#!/bin/bash

#get bedfile for all ORFs that have evidence of translation OR are canonical and longer than 25 nt, save orf information in all.bed
Rscript getbed_all.R

#find which ORFs overlap each other 
bedtools intersect -a /home/aar75/rna_seq/Salmon_20221011/all.bed -b /home/aar75/rna_seq/Salmon_20221011/all.bed -wo -s > overlap_info.txt
#get list of ORFs to remove, specifically find ORF pairs that overlap on the same strand, For ORF pairs where overlap of either ORF is >=0.75, 
#remove the shorter orf (ie get its name) which is stored in /home/aar75/rna_seq/Salmon_20221011/overlapORFs2remove.RDS, 710 orfs
Rscript getOrfs2remove.R

# Rscript getBedFasta.R gets bedfile for annotated and unannotated ORFs and gets transcript sequences for annotated genes (saved as annotated.fasta)
# 25,963 ORFs
Rscript getBedFasta.R

#get transcript sequences for noncanonical orfs
bedtools getfasta -fi /home/aar75/rna_seq/Salmon_20221011/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa -name -s -bed /home/aar75/rna_seq/Salmon_20221011/unannotated.bed -fo temp.fa

#remove the (+) or )-) that bedtools adds to transcript names
awk '{gsub(/\(\+\)|\(\-\)/, "", $0); print}' temp.fa > unannotated.fa 

#concatenate fasta files
cat unannotated.fa annotated.fasta > transcriptome.fa

#get salmon index 
salmon index -t transcriptome.fa -i k_25_idx -k 25 #with kmer size =25

#salmon then removes 31 ORFs bc they are exact sequence duplicates

Rscript getFinalORFlist.R
#so final transcriptome list = #24,570
#also adds transcriptome list to sql table 'expression_transcriptome', 

# conserved    denovo     other
#      5638      8116     10816
#                   Dubious                      None                pseudogene
#                       228                     18455                         5
# transposable_element_gene           Uncharacterized                  Verified
#                        30                       710                      5142
   # canonical noncanonical
   #      5882        18688



#get txt files that list samples to loop through (ie samples_PE.txt and samples_SE.txt)
Rscript getsamplelists.R 

#run alignments
./trim_align_PE.sh
./trim_align_SE.sh


#get tpm, raw counts and effective lengths
Rscript get_TPM_PE.R
Rscript get_TPM_SE.R

#get alignment info
./get_alignment_info_PE.sh
./get_alignment_info_SE.sh

#clean up orfs and samples from raw expression count files, ie remove samples that are not stranded or have less than 1 million reads
#also adds updates transcriptome list to sql table 'expression_transcriptome', 
#saves finalized expression data into /home/aar75/rna_seq/Salmon_20221011/raw_counts.RDS 
#which is a matrix with 24,514 ORFs x 3,916 samples, containing both PE and SE counts
Rscript cleanup_expression_files.R