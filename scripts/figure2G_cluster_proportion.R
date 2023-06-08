library(tidyverse)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
coexpression_orf_list <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")
cluster_assignments <- data.table::fread("data/interim/wgcna_clr_corr/df_12.0_50_rho_spqn.csv")

expanded_network <- data.table::fread("/home/aar75/coexpression/20221110_rho/networkChanges/networks/whole_999.csv", header = T)


cluster_assignments_network_subset <- cluster_assignments[(cluster_assignments$transcript %in% expanded_network$row | cluster_assignments$transcript %in% expanded_network$column) & (cluster_assignments$orf_cluster != 0)]

cluster_counts <- cluster_assignments_network_subset %>%
    group_by(orf_cluster, is_canonical) %>%
    summarise(n = n()) %>%
    group_by(orf_cluster) %>%
    mutate(total = sum(n), proportion = n / total)

noncanonical_cluster_counts <- cluster_counts[cluster_counts$is_canonical == "noncanonical", ]

noncanonical_cluster_counts$breaks <- cut(noncanonical_cluster_counts$proportion, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
col <- circlize::colorRamp2(c(0, 1), c("#7570b3", "#1b9e77"))
ggplot(noncanonical_cluster_counts, aes(x = proportion)) +
    geom_histogram(binwidth = 0.1, aes(fill = breaks), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    theme_bw() +
    scale_fill_manual(
        values = col(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)),
    ) +
    xlab("Proportion of nORFs in a cluster") +
    ylab("Number of clusters") +
    theme(
        legend.position = "none"
    )
ggsave("reports/figures/paper_figures_01062023/norf_cluster_proportion.pdf", width = 4, height = 4)

col <- circlize::colorRamp2(c(0, 1), c("#CC6677", "#89CCED"))
ggplot(noncanonical_cluster_counts, aes(x = proportion)) +
    geom_histogram(binwidth = 0.1, aes(fill = breaks), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    theme_bw() +
    scale_fill_manual(
        values = col(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)),
    ) +
    xlab("Proportion of nORFs in a cluster") +
    ylab("Number of clusters") +
    theme(
        legend.position = "none",
        axis.text=element_text(size=16),
        axis.title = element_text(size=24),
    )
ggsave("reports/figures/paper_figures_01062023/norf_cluster_proportion_talk.png", width = 8, height = 8)
