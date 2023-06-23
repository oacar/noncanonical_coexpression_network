#script to plot figure 3D nORFs and cORFs that are Sfp1 or Hsf1 targets are more likely to be downregulated when Sfp1 or Hsf1 are deleted compared to ORFs that are not targets 
#RNA seq data is from SRA accessions: SRP437124 (Hsf1) and SRP159150 (Sfp1)

library(DESeq2)
library(ggplot2)
library(dplyr)
library(RMariaDB)
library(igraph)
library(ggpubr)
library(EnvStats)
library(rstatix)
#load orf annotation info
#get ORF data

conn <- dbConnect(MariaDB(),groups='mariaDB')
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
gsea<- dbGetQuery(conn, "select * from omer.fgsea_go_slim where padj<0.01;")
dbDisconnect(conn)

orf_info%>%filter(transcript%in%gsea$transcript[gsea$pathway!='GO:0008150'])%>%count(is_canonical)
orf_info%>%filter(transcript%in%colnames(rho))%>%count(is_canonical)
any(gsea$pathway=='GO:0008150')

gsea[gsea$pathway!='GO:0008150' & gsea$transcript%in%orf_info$transcript[orf_info$is_canonical=='noncanonical'],]%>%group_by(transcript)%>%count()%>%arrange((n))



rho<-readRDS('spqn_raw5_sample400.RDS')

#import TF dist matrix
orf_tf_dist<-readRDS('ORF_TF_distances_TSS_all.RDS')
cutoff<-200
orf_tf_dist[orf_tf_dist<=cutoff] <-1
orf_tf_dist[orf_tf_dist>cutoff] <-0
orf_tf_dist[is.na(orf_tf_dist)]<-0

orf_tf_dist<-as.data.frame(orf_tf_dist)

orf_tf_dist$transcript<-rownames(orf_tf_dist)
orf_tf_dist<-filter(orf_tf_dist,transcript %in% rownames(rho))
orf_tf_dist<-left_join(orf_tf_dist,orf_info[,c('transcript','classification','is_canonical','orf_class')])
rownames(orf_tf_dist)<-orf_tf_dist$transcript

# import data for Sfp1
raw_data<-readRDS('/home/aar75/rna_seq/Salmon_20221011/raw_counts.RDS') 
raw_data<-raw_data[rownames(orf_tf_dist),]

#import study list
studylist<-read.csv('SRP159150..txt',stringsAsFactors=FALSE)
alignment_info<-read.delim('alignment_info_5_9_22_PE.txt',header=FALSE,sep = '\t',stringsAsFactors=FALSE,col.names = c('sample','num_reads','map_rate','library_type'))

alignment_info$reads_mapped<-((alignment_info$map_rate/100)*alignment_info$num_reads)
samples2remove<-alignment_info$sample[which(alignment_info$library_type=="IU"|
                                              alignment_info$library_type==""|
                                              alignment_info$library_type=="U"|
                                              alignment_info$reads_mapped <1000000)]



#filter study list to contain only samples from SRP159150

studylist<- filter(studylist, Organism=='Saccharomyces cerevisiae' & Genotype %in% c('wild type', 'delta sfp1') & growth_phase =='log')
studylist$genotype<-gsub(x=studylist$Genotype,pattern = ' ',replacement = '',ignore.case = TRUE)
#studylist <- studylist[order(studylist$sample),]

idx<-which(colnames(raw_data) %in% studylist$Run)
raw<-raw_data
raw<-raw[,idx]

rownames(studylist)<-studylist$Run
studylist$group<-factor(studylist$genotype)

raw<-as.matrix(raw)
raw<-round(raw)
#all(rownames(studylist) == colnames(raw))
dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = studylist,
                              design = ~ group)
dds <- DESeq(dds)

rld <- rlog(dds, blind=TRUE)
res <- results(dds, contrast=c("group","deltasfp1","wildtype"))
res<-as.data.frame(res)


res$transcript<-rownames(res)
res<-left_join(res,orf_tf_dist[,c('transcript','is_canonical')])

TFBS_orfs<-filter(orf_tf_dist, Sfp1 == 1)$transcript
res<-mutate(res, TF_binding=case_when(transcript %in% TFBS_orfs ~ 'yes', 
                                      TRUE ~'no'))


res<-mutate(res, result= case_when(log2FoldChange>0.5 & padj<0.05~ 'upregulated', 
                                   log2FoldChange < (-0.5) & padj <0.05 ~ 'downregulated', 
                                   is.na(padj) ~ 'not expressed', 
                                   TRUE ~ 'not DE'))
res<-mutate(res, result2=case_when(result=='downregulated'~'yes',
                                   TRUE ~'no'))
cORF<-filter(res, is_canonical=='canonical')
nORF<-filter(res, is_canonical=='noncanonical')

cm_c<-table(cORF$result2, cORF$TF_binding)
cm_n<-table(nORF$result2, nORF$TF_binding)
sfp1_counts_canonical<-cm_c
sfp1_counts_noncanonical<-cm_n


###### import rna seq raw counts for Hsf1
#import study list
studylist<-read.csv('SRP437124.csv',stringsAsFactors=FALSE)
alignment_info<-read.delim('/home/aar75/rna_seq/SRP437124/alignment_info_PE.txt',header=FALSE,sep = '\t',stringsAsFactors=FALSE,col.names = c('sample','num_reads','map_rate','library_type'))

alignment_info$reads_mapped<-((alignment_info$map_rate/100)*alignment_info$num_reads)
samples2remove<-alignment_info$sample[which(alignment_info$library_type=="IU"|
                                              alignment_info$library_type==""|
                                              alignment_info$library_type=="U"|
                                              alignment_info$reads_mapped <1000000)]

raw<-read.csv('/home/aar75/rna_seq/SRP437124/RAW_counts_PE.csv',stringsAsFactors = FALSE)
rownames(raw)<-raw$Name
raw$Name<-NULL
raw<-raw[rownames(orf_tf_dist),]


studylist<-filter(studylist, condition %in% c("hsf1del_hs","wt_hs"))
samples<-intersect(colnames(raw),studylist$sample)
studylist<-filter(studylist, sample %in% samples)
raw<-raw[,studylist$sample]

#do differential expression:
rownames(studylist)<-studylist$sample
studylist$group<-factor(studylist$condition)

raw<-as.matrix(raw)
raw<-round(raw)
all(rownames(studylist) == colnames(raw))
dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = studylist,
                              design = ~ group)
dds <- DESeq(dds)

rld <- rlog(dds, blind=TRUE)




res <- results(dds, contrast=c("group","hsf1del_hs","wt_hs"))

res<-as.data.frame(res)

res$transcript<-rownames(res)
res<-left_join(res,orf_tf_dist[,c('transcript','is_canonical')])

TFBS_orfs<-filter(orf_tf_dist, Hsf1 == 1)$transcript
res<-mutate(res, TF_binding=case_when(transcript %in% TFBS_orfs ~ 'yes', 
                                      TRUE ~'no'))


res<-mutate(res, result= case_when(log2FoldChange>0.5 & padj<0.05~ 'upregulated', 
                                   log2FoldChange < (-0.5) & padj <0.05 ~ 'downregulated', 
                                   is.na(padj) ~ 'not expressed', 
                                   TRUE ~ 'not DE'))

res<-mutate(res, result2=case_when(result=='downregulated'~'yes',
                                   TRUE ~'no'))


cORF<-filter(res, is_canonical=='canonical')
nORF<-filter(res, is_canonical=='noncanonical')


#association between an ORF having a TF binding event (in chip Exo) and being downregulated when TF is deleted (compared to WT)
cm_c<-table(cORF$result2, cORF$TF_binding)
cm_n<-table(nORF$result2, nORF$TF_binding)

hsf1_counts_canonical<-cm_c
hsf1_counts_noncanonical<-cm_n

odds_ratio_df <- data.table::rbindlist(
  list(  pairwise_fisher_test(hsf1_counts_canonical,detailed = T)%>%mutate(tf ='hsf1',cat='canonical'),
         pairwise_fisher_test(hsf1_counts_noncanonical,detailed = T)%>%mutate(tf ='hsf1',cat='noncanonical'),
         pairwise_fisher_test(sfp1_counts_canonical,detailed = T)%>%mutate(tf ='sfp1',cat='canonical'),
         pairwise_fisher_test(sfp1_counts_noncanonical,detailed = T)%>%mutate(tf ='sfp1',cat='noncanonical')
  )) 


significance_size <- 4 
axis_title_size = 12
legend_text_size = 10
legend_title_size = 10
axis_text_size <- 10

odds_ratio_df%>%
  ggplot(aes(x=tf, y=estimate,fill=cat))+
  geom_col(position=position_dodge2())+
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width=.2, color='black',size=.5,
                position= position_dodge(width=.9))+
  theme_bw()+
  xlab("Transcription factor")+
  ylab("Likelihood of down regulation\nupon TF knock-out\nTF targets vs non-targets\n(odds ratio)")+
  scale_fill_manual(name='',labels=c('canonical','noncanonical'), values = c('#7772b4',"#219f78" ))+
  scale_x_discrete(labels=c(expression(Delta*'Hsf1'), expression(Delta*'Sfp1')))+
  annotate("text", x = 1, y = 34, label = "ns", size = 4, color = "black") +
  annotate("text", x = 2, y = 34, label = "ns", size = 4, color = "black") +
  theme(legend.position='bottom')+geom_hline(yintercept=1, linetype='dashed')+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) +
  font("legend.text", size=legend_text_size)+
  font("legend.title", size=legend_title_size)
ggsave('~/coexpression/tf_deletion_odds.pdf',height = 3, width=4)


