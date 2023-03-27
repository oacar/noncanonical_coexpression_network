#!/usr/bin/Rscript
#clean up orfs and samples from expression files ie remove unstranded rna seq samples and remove samples with less than 1 million reads
library(dplyr)
library(RMariaDB)

#read in expression data
raw_SE<-read.csv('RAW_counts_SE.csv')
raw_PE<-read.csv('RAW_counts_PE.csv')

#remove ORFs from expression files
raw<-left_join(raw_PE,raw_SE)
rm(raw_SE,raw_PE)

#read in alignment info
alignment_info_SE<-read.delim('alignment_info_SE.txt',header=FALSE,sep = '\t',col.names = c('sample','num_reads','map_rate','library_type'))
alignment_info_PE<-read.delim('alignment_info_PE.txt',header=FALSE,sep = '\t',col.names = c('sample','num_reads','map_rate','library_type'))

alignment_info<-rbind(alignment_info_PE,alignment_info_SE)
rm(alignment_info_PE,alignment_info_SE)

alignment_info$reads_mapped<-((alignment_info$map_rate/100)*alignment_info$num_reads)
length(unique(alignment_info$sample)) #4005
alignment_info<-filter(alignment_info, library_type %in% c('SR','SF','ISR','ISF'))
length(unique(alignment_info$sample)) #3925
alignment_info<-filter(alignment_info, reads_mapped >= 1e+06) 
length(unique(alignment_info$sample)) #3916

rownames(raw)<-raw$Name
raw$Name<-NULL

#remove samples
idx<-which(colnames(raw) %in% alignment_info$sample)
ncol(raw) #4005
raw<-raw[,idx]
ncol(raw) #3916

dim(raw)
#24514  3916
saveRDS(raw,'raw_counts.RDS')



