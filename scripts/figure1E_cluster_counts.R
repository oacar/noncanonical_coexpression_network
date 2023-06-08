library(tidyverse)
con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")

cluster_df <- data.table::fread("/home/oma21/coexpression/data/interim/wgcna_clr_corr/df_12.0_50_rho_spqn.csv")

cluster_go_counts <- data.table::fread("/home/oma21/coexpression/data/interim/wgcna_clr_corr//cluster_group_counts_all_12.0_50_rho_spqn.csv")

cluster_go_terms <- data.table::fread("/home/oma21/coexpression/data/interim/wgcna_clr_corr//cluster_group_counts_all_12.0_50_rho_spqn_go.csv")
cluster_go_counts <- cluster_go_counts %>% left_join(
    cluster_df %>%
        group_by(orf_cluster) %>%
        count(is_canonical) %>%
        ungroup() %>%
        pivot_wider(names_from = is_canonical, values_from = n) %>%
        mutate(perc_canonical = canonical / (canonical + noncanonical))
)
cluster_go_counts$go_enriched <- ifelse(cluster_go_counts$bp_count > 0, "GO enriched", "Not GO enriched")

ggplot(
    cluster_go_counts %>%
        filter(canonical > 5) %>%
        group_by(go_enriched) %>%
        summarise(n = n()),
    aes(x = go_enriched, y = n)
) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.5, hjust = 0.5, size = 3) +
    theme_classic() +
    labs(x = "", y = "Number of clusters")
ggsave("reports/figures/paper_figures_01062023/cluster_go_canonical.pdf", width = 3, height = 3)
