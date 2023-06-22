#script to plot supp figure 13B: boxplots of orientation vs biological process similarity
#omer.gsea_similarity_random data is generated in gsea_jaccard_random.r
# omer.fgseas_slim_similarity data is generated in gsea_relevance_similarity.R

library(dplyr)
library(magrittr)
library(stringr)
library(ggpubr)
library(forcats)
library(scales)
library(EnvStats)
library(RMariaDB)
library(effsize)

#path to save plots
file_path<-'./20230106_figures/'


orientation_levels <- c("background", "down\nopposite", "down\nsame", "antisense\noverlap", "up\nopposite", "up\nsame")
significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


conn <- dbConnect(MariaDB(), 
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
background_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_bp_noexponent_similarity_random')
orf_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_similarity')
dbDisconnect(conn)

background_pairs<-data.frame('orientation'='background', 'rel'=background_pairs)

orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down same',
                                                   orientation=='upsame' ~'up same',
                                                   orientation =='antisense-overlap' ~'antisense overlap',
                                                   orientation=='downopposite' ~'down opposite',
                                                   orientation=='upopposite' ~'up opposite'))

#import coexpression matrix
rho_matrix<-readRDS('spqn_raw5_sample400.RDS')
numObs<-readRDS('numobs_raw5_sample400.RDS')
rho_matrix[numObs<400]<-NA


orf_pairs$rho<-NA
orf_pairs<-filter(orf_pairs, distance <=500)
orf_pairs<-filter(orf_pairs, primary_orf %in% rownames(rho_matrix))
orf_pairs<-filter(orf_pairs, neighbor_orf %in% rownames(rho_matrix))
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho_matrix[orf_pairs$neighbor_orf[i],orf_pairs$primary_orf[i]]
}
idx<- which(is.na(orf_pairs$rho) & !is.na(orf_pairs$rel))
orf_pairs$rel[idx]<-NA
idx<-which(is.na(orf_pairs$rel) & !is.na(orf_pairs$rho))
orf_pairs$rel[idx]<-0
orf_pairs<-orf_pairs[,c('orientation','rel')]

orf_pairs<-rbind(orf_pairs,background_pairs)
orf_pairs<-filter(orf_pairs, !is.na(rel))

pdf(sprintf("%sfigure2_orientationVsBPsimilarity.pdf", file_path), width = 3.6, height = 3)
orf_pairs %>%
  mutate(orientation = factor(str_replace(orientation, "\ ", "\n"), levels = orientation_levels)) %>%
  ggboxplot(x = "orientation", y = "rel", ylab = "Biological process enrichment similarity", notch = TRUE, fill = "#E59F01") +
  geom_hline(yintercept=median(filter(orf_pairs,orientation=='background')$rel),linetype='dotted')+
  stat_n_text(size = 2.75) + 
  stat_compare_means(label = "p.signif", ref.group = "background", label.y = 1.1, size = significance_size) +
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  rremove("xlab")
dev.off()


#down same
o<-'down same'
wilcox.test(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
cliff.delta(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
#down same: large (Cliff’s Delta d=0.501, Mann-Whitney U-test, p< 2.2e-16)

o<-'up same'
wilcox.test(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
cliff.delta(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
#up same: small (Cliff’s Delta d=0.169, Mann-Whitney U-test, p< 2.2e-16)

o<-'antisense overlap'
wilcox.test(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
cliff.delta(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
#antisense overlap: negligible (Cliff’s Delta d=-0.0359, Mann-Whitney U-test, p=0.163)

o<-'up opposite'
wilcox.test(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
cliff.delta(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
#up opposite: negligible (Cliff’s Delta d=0.0573, Mann-Whitney U-test, p=0.0169)

o<-'down opposite'
wilcox.test(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
cliff.delta(filter(orf_pairs,orientation==o)$rel, filter(orf_pairs,orientation=='background')$rel)
#down opposite: negligible (Cliff’s Delta d=-0.015, Mann-Whitney U-test, p=0.469)

