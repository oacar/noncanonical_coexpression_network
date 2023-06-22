#this script plots figure 3C: looking at the association between being coexpressed with genes involved in transport and having a predicted transmembrane domain
#x axis is whether an nORF has a significant GSEA enrichment for any 15 GO slim transport terms: transport, ion transport, a
#amino acid transport, lipid transport, carbohydrate transport, regulation of transport, transmembrane transport, vacuolar transport, vesicle-mediated transport,
#endosomal transport, nucleobase-containing compound transport, Golgi vesicle transport, nucleocytoplasmic transport, nuclear transport, or 
#cytoskeleton-dependent intracellular transport.

#y axis is number of nORFs that have a predicted transmembrane (TM) domains. TM domains were predicted using TMHMM 2.0

#gsea enrichments (ie fgsea_slim_bp_noexponent.csv) for each ORF generated in file run_gsea.R 

library(dplyr)
library(RMariaDB)
library(ggplot2)

gseas<-read.csv('fgsea_slim_bp_noexponent.csv')

gos<-unique(gseas$pathway)


conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
TM_info<-dbGetQuery(conn, "SELECT * FROM carly.carly_orfs_TMs")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_ids<-dbGetQuery(conn, "SELECT * FROM CARVUNIS_YEAST.orf_ids_crossmap")
dbDisconnect(conn)

orf_ids$transcript<-sprintf('chr%i_%i',orf_ids$chr_num,orf_ids$coor1)

#get list of genes involved in transport
go_list <- read.delim("go_slim_mapping.tab", stringsAsFactors = F, header = F)
colnames(go_list) <- c("gene", "common_name", "sgd_id", "go_type", "goTerm", "go_id", "orf_type")


#get list of nORFs enriched for transport 
go_list<-unique(go_list[,c('goTerm','go_id')])
x<-left_join(gseas,go_list, by=c('pathway'='go_id'))

rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

# denovo_df<-filter(orf_info, classification=='denovo' & transcript %in% x$transcript)
nORF_df<-filter(orf_info, is_canonical=='noncanonical' & transcript %in% x$transcript & transcript %in% rownames(rho))
nORF_df<-left_join(nORF_df, TM_info)
nORF_df<-mutate(nORF_df, TM=case_when(Number_of_predicted_TMHs ==0 ~ 'no',
                                      Number_of_predicted_TMHs>0 ~'yes'))


nORF_df$transport<-NA
for(i in 1:nrow(nORF_df)){
  y<-filter(x, transcript==nORF_df$transcript[i] & padj<0.01)
  if(nrow(y)>0){
    idx<-grep(x=y$goTerm, pattern='transport')
    if(length(idx)>0){nORF_df$transport[i]<-'yes'}
    else{nORF_df$transport[i]<-'no'}
  }
}

cm_nORF<-table(nORF_df$TM, nORF_df$transport)
fisher.test(cm_nORF)


library(rstatix)
library(tibble)
library(ggpubr)
rho_hex<-'#DC0000B2'
file_path<-'./20230106_figures/'
significance_size <- 4 
axis_title_size = 12
legend_text_size = 10
legend_title_size = 10
axis_text_size <- 10


getPlotdf<-function(cm){
  tp = 100*(cm[2,2]/(cm[2,2]+cm[1,2]))
  ntp = 100*(cm[2,1]/(cm[1,1]+cm[2,1]))
  plot_df<-data.frame('transport'=c('transport','not transport'),'percentage'=c(tp, ntp), 'denom'=c((cm[2,2]+cm[1,2]),(cm[1,1]+cm[2,1])))
  return(plot_df)
}


complex_df<-getPlotdf(cm_nORF)
complex_stats<-pairwise_fisher_test(cm_nORF) %>% add_column('transport'='yes')

complex_df <- complex_df %>%
  mutate(se= sqrt(percentage*(100-percentage)/(denom)))

pdf(sprintf('%sTM_enrichment_nORF.pdf',file_path),width=2, height = 3)
complex_df %>% mutate(transport =case_when(transport=='transport' ~'yes',
                                           transport=='not transport' ~'no')) %>%
  ggbarplot( x="transport", y="percentage", xlab="coexpressed\nwith transport\nrelated genes",ylab ="% of nORFs with a predicted TM domain",
             fill = "transport", color=NA, palette = c('#BFBADA','#E32726'),
             label = FALSE)+
  geom_errorbar(aes(group=transport, ymax = percentage + se, ymin=percentage - se),
                position = position_dodge(width = 0.8), width = 0.25)+
  stat_pvalue_manual(complex_stats, x="transport", label="p.adj.signif",y.position=21, size = significance_size)+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + grids(axis = c("xy", "x", "y"), color = "grey92", size = NULL, linetype = NULL)+
  rremove("legend")
dev.off()
