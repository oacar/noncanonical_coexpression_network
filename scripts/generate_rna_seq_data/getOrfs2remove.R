library(dplyr)

#for ORFs that overlap on the same strand - keep the longer ORF

overlap<-read.delim('overlap_info.txt',col.names=c('orf1_chr','orf1_coord1','orf1_coord2','orf1','orf1_x','orf1_strand',
                                                                                       'orf2_chr','orf2_coord1','orf2_coord2','orf2','orf2_x','orf2_strand','overlap_amount'))
overlap$orf1_length<-overlap$orf1_coord2-overlap$orf1_coord1
overlap$orf2_length<-overlap$orf2_coord2-overlap$orf2_coord1
overlap$orf1_overlap_prop<-overlap$overlap_amount/overlap$orf1_length
overlap$orf2_overlap_prop<-overlap$overlap_amount/overlap$orf2_length

overlap<-overlap[,c('orf1','orf2','overlap_amount','orf1_overlap_prop','orf2_overlap_prop','orf1_length','orf2_length')]

overlap<-filter(overlap, orf1 != orf2)

highoverlap<-filter(overlap,orf1_overlap_prop>=0.75 | orf2_overlap_prop>=0.75)
nrow(highoverlap) #2114


orfs2remove<-unique(c(filter(highoverlap, orf1_length>orf2_length)$orf2,filter(highoverlap,orf1_length<orf2_length)$orf1))

saveRDS(orfs2remove,'overlapORFs2remove.RDS') #710 orfs to remove 
