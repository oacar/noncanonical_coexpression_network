#!/usr/bin/env Rscript

#combine all the rho or num obs columns into one matrix
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
input_files_path<-args[1]
output_file_name<-args[2]

files<-list.files(path = input_files_path, full.names = TRUE, recursive = FALSE)
x<-read.delim(files[1],stringsAsFactors = FALSE,sep = ' ')
for(i in 2:length(files)){
  temp<-read.delim(files[i],stringsAsFactors = FALSE,sep = ' ')
  x<-left_join(x,temp)
}

rownames(x)<-x$transcript
x$transcript<-NULL
saveRDS(x,output_file_name)