#plots network transitivity for canonical network, expanded network and 1000 randomized networks
#1000 randomized networks generated in script rewire_network.py
#network transitivity for 1000 randomized networks calculated in get_network_properties.R
#canonical and expanded edge list generated in get_networks_for_plotting.R


library(dplyr)
library(igraph)
library(ggplot2)

significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


tt<-read.delim('transitivity_random.txt',header=FALSE)
plot_df<-data.frame('transitivity'=tt[,1],'network_type'='random')

#import real networks
expanded<-read.csv('expanded_network.csv')
canonical<-read.csv('canonical_network.csv')


expanded_network<-graph_from_data_frame(expanded[,c('row','column')],directed = FALSE) 
temp<-data.frame('transitivity'=transitivity(expanded_network),
                 'network_type'='expanded')
plot_df<-rbind(plot_df,temp)


canonical_network<-graph_from_data_frame(canonical[,c('row','column')],directed = FALSE) 
temp<-data.frame('transitivity'=transitivity(canonical_network),
                 'network_type'='canonical')
plot_df<-rbind(plot_df,temp)
plot_df$network_type<-factor(plot_df$network_type, levels = c('canonical','expanded','random'))



pdf("./20230106_figures/transitivityPlot.pdf", width = 4.5, height = 2)  
ggplot() +
  geom_density(data = plot_df %>% filter(network_type== "random"), aes(x = transitivity, y = ..scaled.., color = network_type),linewidth=0.6) +
  geom_vline(data =plot_df %>% filter(network_type == "canonical"), aes(xintercept = transitivity, color = network_type),linewidth=0.6) +
  geom_vline(data =plot_df %>% filter(network_type == "expanded"), aes(xintercept = transitivity, color = network_type),linewidth=0.6) +
  guides(color = guide_legend(title = "Network type"), shape = guide_legend(title = "Network type")) +
  ylab("Density") +
  xlab("network transitivity") +
  theme_bw() +
  theme(panel.spacing = unit(2, "lines"), legend.position = "bottom",axis.text=element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size), legend.text = element_text(size=legend_text_size),
        legend.title = element_text(size=legend_text_size),legend.key.size=unit(0.1, "lines")) +
  scale_color_manual(values = c('#7671B2', '#F26322','#FCAF17')) +
  scale_y_continuous(expand = c(0, 0))
dev.off()
