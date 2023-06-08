library(tidyverse)
library(igraph)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")

coexpression_orf_list <- DBI::dbGetQuery(con, "select * from coexpressionOrfList_blaste4")

mat_matrix <- readRDS("/home/aar75/coexpression/20221110_rho/spqn_raw5_sample400.RDS")
obs_matrix <- readRDS("/home/aar75/coexpression/20221110_rho/numobs_raw5_sample400.RDS")
# mat_file <- "/home/oma21/coexpression/data/interim/coexpression_matrix_rho_spqn.csv"
# obs_file <- "/home/oma21/coexpression/data/interim/pairwise_obs.csv"

# mat <- data.table::fread(mat_file)
# obs <- data.table::fread(obs_file)

# mat_matrix <- as.matrix(mat[, 2:ncol(mat)])
# obs_matrix <- as.matrix(obs[, 2:ncol(mat)])

mat_matrix[obs_matrix < 400] <- NA
rownames(obs_matrix) <- colnames(obs_matrix)
rownames(mat_matrix) <- colnames(mat_matrix)
dim(mat_matrix)
matrix_long <- reshape2::melt(mat_matrix) %>% filter(Var1 != Var2)
n_v_c <- c()
n_v_nc <- c()
n_v <- c()
n_e <- c()
quantile_levels <- c(
    .99, .9905, .991, .9915, .992, .9925, .993, .9935, .994, .9945, .995, .9955, .996, .9965, .997, .9975, .998, .9985, .999
)
thresholds <- c()
for (i in quantile_levels) {
    thr <- quantile(mat_matrix[upper.tri(mat_matrix)], i, na.rm = T)
    print(thr)
    thresholds <- c(thresholds, thr)
    g <- matrix_long %>%
        filter(value > thr) %>%
        graph_from_data_frame(directed = F) %>%
        simplify()
    n_v_c <- c(n_v_c, length(V(g)$name[V(g)$name %in% coexpression_orf_list$transcript[coexpression_orf_list$is_canonical == "canonical"]]))
    n_v_nc <- c(n_v_nc, length(V(g)$name[V(g)$name %in% coexpression_orf_list$transcript[coexpression_orf_list$is_canonical == "noncanonical"]]))
    n_v <- c(n_v, length(V(g)))
    n_e <- c(n_e, length(E(g)))
}
res <- data.frame(
    quantile = quantile_levels,
    thresholds = thresholds,
    canonical_count = n_v_c,
    canonical_proportion = n_v_c / length(coexpression_orf_list$transcript[coexpression_orf_list$transcript %in% rownames(mat_matrix) & coexpression_orf_list$is_canonical == "canonical"]),
    noncanonical_count = n_v_nc,
    noncanonical_proportion = n_v_nc / length(coexpression_orf_list$transcript[coexpression_orf_list$transcript %in% rownames(mat_matrix) & coexpression_orf_list$is_canonical == "noncanonical"]),
    total_count = n_v,
    edge_count = n_e
)
p1 <- res %>%
    pivot_longer(cols = c(canonical_proportion, noncanonical_proportion)) %>%
    ggplot(aes(x = quantile, y = value, color = name)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    geom_vline(xintercept = .998, linetype = "dashed") +
    labs(x = "Quantile Threshold", y = "Proportion of ORFs in the network\nafter thresholding") +
    scale_color_manual(
        name = "ORF Type", labels = c("Canonical", "Noncanonical"),
        values = c("#7570b3", "#1b9e77")
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
leg <- cowplot::get_legend(p1)
p2 <- res %>%
    ggplot(aes(x = quantile, y = edge_count)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = .998, linetype = "dashed") +
    labs(x = "Quantile Threshold", y = "Number of Edges in the network\nafter thresholding") +
    theme_bw() +
    scale_y_continuous(breaks = seq(0, 1e6, by = 1e5))
cowplot::plot_grid(p1 + theme(legend.position = "none"), p2, leg, ncol = 2, rel_widths = c(1, 1), rel_heights = c(1, .1))
ggsave("reports/figures/paper_figures_01062023/quantile_threshold_changes.png", width = 9, height = 4, units = "in", dpi = 300, bg = "white")
