#!/usr/bin/env Rscript
library(dplyr)
samples<-read.csv('/home/aar75/rna_seq/Salmon_20221011/9_27_22_studylist.csv')
table(samples$layout)

paired<-filter(samples, layout == 'paired')
single<-filter(samples, layout == 'single')

write.table(paired$sample,file = 'samples_PE.txt',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(single$sample,file = 'samples_SE.txt',quote = FALSE,row.names = FALSE,col.names = FALSE)
