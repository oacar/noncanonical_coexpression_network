#script to plot supplementary figure 5: TF enrichment at various coexpression cutoffs (ie is there an association between an 
#ORF being coexpressed and having their promoters bound by a common TF) using different coexpression cutoffs. 


#contingency matrices to generate this plot is from: get_TF_enrichment_clr_rho_spqn.R

library(ggplot2)
library(dplyr)
library(scales)
library(rstatix)
library(tibble)
library(ggpubr)
rho_hex<-'#DC0000B2'
file_path<-'./20230106_figures/'
significance_size <- 4 
axis_title_size = 12
legend_text_size = 10
legend_title_size = 12
axis_text_size <- 10


getPlotdf_csv<-function(cm){
  cm$X<-NULL
  colnames(cm)<-c('0','1')
  rownames(cm)<-c('0','1')
  cp = (cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = (cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed pairs','not coexpressed pairs'),'prop'=c(cp, ncp),'count'=c(cm['1','1'],cm['1','0'])) plot_df$se<-sqrt(plot_df$prop * (1 - plot_df$prop) / plot_df$count)
  plot_df$prop<-plot_df$prop*100
  plot_df$se<-plot_df$se*100
  return(plot_df)
}
getPlotdf_rds<-function(cm){
  cp = (cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = (cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed','not coexpressed'),'prop'=c(cp, ncp),'count'=c(cm['1','1'],cm['1','0']))
  plot_df$se<-sqrt(plot_df$prop * (1 - plot_df$prop) / plot_df$count)
  plot_df$prop<-plot_df$prop*100
  plot_df$se<-plot_df$se*100
  return(plot_df)
}



make_combined_plot_df<-function(rho_cutoff){
  
  ccm<-readRDS(sprintf('rho_cc_pairswithTFBS_200_%i.RDS',rho_cutoff))
  cc<-getPlotdf_rds(ccm)
  cc$pair_type<-'canonical\ncanonical'
  
  cnm<-readRDS(sprintf('rho_cn_pairswithTFBS_200_%i.RDS',rho_cutoff))
  cn<-getPlotdf_rds(cnm)
  cn$pair_type<-'canonical\nnoncanonical'
  
  nnm<-readRDS(sprintf('rho_nn_pairswithTFBS_200_%i.RDS',rho_cutoff))
  nn<-getPlotdf_rds(nnm)
  nn$pair_type<-'noncanonical\nnoncanonical'
  
  plot_df<-rbind(cc,cn)
  plot_df<-rbind(plot_df,nn)
  plot_df$coexp_cutoff<-rho_cutoff
  return(plot_df)
}

make_stats_df<-function(rho_cutoff){
  
  ccm<-readRDS(sprintf('rho_cc_pairswithTFBS_200_%i.RDS',rho_cutoff))
  cnm<-readRDS(sprintf('rho_cn_pairswithTFBS_200_%i.RDS',rho_cutoff))
  nnm<-readRDS(sprintf('rho_nn_pairswithTFBS_200_%i.RDS',rho_cutoff))
  
  stats_cc<-pairwise_fisher_test(ccm, p.adjust.method = "holm",detailed=TRUE) %>% add_column('pair_type'='canonical\ncanonical')
  #cc odds ratio = 4.2848 p-value <2.2e-16
  stats_cn<-pairwise_fisher_test(cnm, p.adjust.method = "holm") %>% add_column('pair_type'='canonical\nnoncanonical')
  #cn odds ratio = 2.9981 p-value < 2.2e-16
  stats_nn<-pairwise_fisher_test(nnm, p.adjust.method = "holm") %>% add_column('pair_type'='noncanonical\nnoncanonical')
  #nn odds ratio = 3.8575 p-value <2.2e-16
  stats_tf<- bind_rows(stats_cc,stats_cn,stats_nn)
  stats_tf<- stats_tf %>% add_column('coexp_cutoff'=rho_cutoff)
  return(stats_tf)
}

plots_df<-make_combined_plot_df(rho_cutoff=90)
plots_df<-rbind(plots_df, make_combined_plot_df(rho_cutoff=95))
plots_df<-rbind(plots_df, make_combined_plot_df(rho_cutoff=99))
plots_df<-rbind(plots_df, make_combined_plot_df(rho_cutoff=999))
stats_df<-make_stats_df(rho_cutoff=90)
stats_df<-rbind(stats_df,make_stats_df(rho_cutoff=95))
stats_df<-rbind(stats_df,make_stats_df(rho_cutoff=99))
stats_df<-rbind(stats_df,make_stats_df(rho_cutoff=999))



p<-plots_df %>% #mutate(coexpressed =case_when(coexpressed=='coexpressed pairs' ~'coexpressed',
  # coexpressed=='not coexpressed pairs' ~'not\ncoexpressed')) %>%
  ggbarplot( x="pair_type", y="prop", xlab="ORF pair type",ylab ="% of pairs with their promoter bound by a common TF",
             fill = "coexpressed", color=NA, palette = c('#3FC1C9','#E32726'),
             label = plots_df$count,   legend = "bottom",
             position = position_dodge(0.8)) + 
  # label = plots_df$count) + 
  # geom_text(aes(y = ..count.., group = interaction(coexpressed, coexp_cutoff)), 
  #           size = 2.5, vjust = -0.5, hjust = 0.5, position = position_dodge(width = .9)) +
  # geom_text(aes(label = format(plot_df_stats$count, big.mark = ",")), size = 2.5, vjust = -0.5,hjust=0.5) +
  stat_pvalue_manual(stats_df, x='pair_type', label="p.adj.signif",y.position=42, size = significance_size)+
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "black") +
  font("x.text", size = 9) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + grids(axis = c("xy", "x", "y"), color = "grey92", size = NULL, linetype = NULL)
pdf(sprintf('%sTFBS_Enrichment_supplementary_figure.pdf',file_path),width=6.75, height = 5)
facet(p,facet.by="coexp_cutoff" ,nrow = 2)
dev.off()





