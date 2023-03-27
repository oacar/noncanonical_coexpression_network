library(RMariaDB)
library(dplyr)
library("Biostrings")
library(stringr)

#get orf and gene list
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
coexpression_orf_list<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


coexpression_orf_list$length<-coexpression_orf_list$coor2-coexpression_orf_list$coor1

#filter to keep orfs that have evidence of translation OR are canonical 
coexpression_orf_list<-filter(coexpression_orf_list, orf_class %in% c('Verified','Uncharacterized','transposable_element_gene') | qval<0.05)
coexpression_orf_list<-filter(coexpression_orf_list, length>=25)

#remove orfs with overlap
orfs2remove<-readRDS('overlapORFs2remove.RDS')
coexpression_orf_list<-filter(coexpression_orf_list, !(transcript %in% orfs2remove)) #24623


annotated<-filter(coexpression_orf_list,orf_class!='None') #6,168
unannotated<-filter(coexpression_orf_list, orf_class=='None') #18,455

#bed files 
#replace chr # with NCBI chromosome IDs
NCBI_chrID<-c('ref|NC_001133|','ref|NC_001134|','ref|NC_001135|','ref|NC_001136|','ref|NC_001137|','ref|NC_001138|',
              'ref|NC_001139|','ref|NC_001140|','ref|NC_001141|','ref|NC_001142|','ref|NC_001143|','ref|NC_001144|',
              'ref|NC_001145|','ref|NC_001146|','ref|NC_001147|','ref|NC_001148|')
#create unannotated bed
unannotated$chr<-NA
for(i in 1:length(NCBI_chrID))
{
  idx<-which(unannotated$chr_num==i)
  unannotated$chr[idx]<-NCBI_chrID[i]
}


idx<-which(unannotated$strand=='+')
unannotated$coor1[idx]<-unannotated$coor1[idx]-1

#make unannotated bed file:
bed_unannotated<-data.frame('chr'=unannotated$chr,'start'=unannotated$coor1,
                            'end'=unannotated$coor2,'name'=unannotated$transcript,'dot'=1,
                            'strand'=unannotated$strand, stringsAsFactors = FALSE)

write.table(bed_unannotated,file = 'unannotated.bed',quote = FALSE,row.names = FALSE,sep="\t",col.names = FALSE)

#make annotated bed file
annotated$chr<-NA
for(i in 1:length(NCBI_chrID))
{
  idx<-which(annotated$chr_num==i)
  annotated$chr[idx]<-NCBI_chrID[i]
}


idx<-which(annotated$strand=='+')
annotated$coor1[idx]<-annotated$coor1[idx]-1
bed_annotated<-data.frame('chr'=annotated$chr,'start'=annotated$coor1,
                          'end'=annotated$coor2,'name'=annotated$transcript,'dot'=1,
                          'strand'=annotated$strand, stringsAsFactors = FALSE)
write.table(bed_annotated,file = 'annotated.bed',quote = FALSE,row.names = FALSE,sep="\t",col.names = FALSE)

##create annotated transcript sequence file
fastaFile <- readAAStringSet("orf_coding_all_R64-2-1_20150113.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
df$seq_name<-str_extract(df$seq_name, pattern='[[:alnum:]-]*')

df<-df[which(df$seq_name %in% annotated$gene ),]
df<-left_join(df,annotated[,c('gene','transcript')],by=c('seq_name'='gene'))

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"transcript"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

writeFasta(df,'annotated.fasta')

