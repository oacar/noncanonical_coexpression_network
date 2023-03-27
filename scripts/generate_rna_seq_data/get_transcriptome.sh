#!/bin/bash

#get bedfile for all ORFs that have evidence of translation OR are annotated, save orf information in all.bed
Rscript getbed_all.R

#find which ORFs overlap each other 
bedtools intersect -a all.bed -b all.bed -wo -s > overlap_info.txt
#get list of ORFs to remove, specifically find ORF pairs that overlap on the same strand, For ORF pairs where overlap of either ORF is >=0.75, 
#remove the shorter orf (ie get its name) which is stored in overlapORFs2remove.RDS
Rscript getOrfs2remove.R

# Rscript getBedFasta.R gets bedfile for unannotated and annotated ORFs and gets transcript sequences for annotated genes (saved as annotated.fasta)
# also adds transcriptome list to sql table 'coexpression_transcriptome'

Rscript getBedFasta.R

#get transcript sequences for unannotated orfs
bedtools getfasta -fi S288C_reference_sequence_R64-2-1_20150113.fsa -name -s -bed unannotated.bed -fo temp.fa

#remove the (+) or )-) that bedtools adds to transcript names
awk '{gsub(/\(\+\)|\(\-\)/, "", $0); print}' temp.fa > unannotated.fa 

#concatenate fasta files
cat unannotated.fa annotated.fasta > transcriptome.fa

