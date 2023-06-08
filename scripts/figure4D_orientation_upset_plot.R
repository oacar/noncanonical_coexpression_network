library(ComplexUpset)
library(tidyverse)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")
tpms <- DBI::dbGetQuery(con, "select * from april.coexpression_median_tpm")
neighbor_info_denovo = DBI::dbGetQuery(con, "select * from omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")

neighbor_info_denovo <- neighbor_info_denovo %>%
    left_join(
        allORFs[, c("transcript", "coor1", "strand")],
        by = c("primary_orf" = "transcript")
    ) %>%
    left_join(
        allORFs[, c("transcript", "coor1")],
        by = c("neighbor_orf" = "transcript")
    ) %>%
    select(primary_orf, neighbor_orf, distance, orientation)

orientation_overlaps <- neighbor_info_denovo %>%
    # filter(classification=='denovo')%>%
    mutate(exists = ifelse(distance < 500, T, F)) %>%
    pivot_wider(names_from = orientation, id_cols = c(primary_orf), values_fill = F, values_from = exists, values_fn = any) %>%
    # rename(`anti-sense overlap` = overlap)%>%
    mutate(independent = (!(`upsame` | `upopposite` | `downsame` | `downopposite` | `antisense-overlap`)) & (primary_orf %in% far_away$transcript))


upset(orientation_overlaps %>% filter(primary_orf %in% tpms$transcript, primary_orf != "chr5_473485"),
    c("downsame", "upsame", "downopposite", "upopposite", "independent", "antisense-overlap"),
    sort_intersections = "ascending",
    base_annotations = list(
        "Intersection size" = intersection_size(

            #  text_mapping =aes(fill=expressed)
        )
    ),
    sort_intersections_by = c("degree", "cardinality"), max_degree = 2, min_degree = 1, min_size = 20, height_ratio = .6, width_ratio = .1, name = "orientation"
)
ggsave("../reports/figures/paper_figures_01062023/orientation_upset_plot.pdf", width = 9, height = 4)
