cn<-colnames(orf_info)


orf_info$source <- ifelse(orf_info$orf_class=='None',"Wacholder2023","SGD")
attribute_df <-orf_info
attribute_df<-attribute_df%>%select(ID=orf_coord_id,orf_id, evolutionary_class = classification, is_canonical, systematic_gene_name = gene, NAME=GENENAME, SGD_orf_classification = orf_class, wacholder2023_qval=qval)
attribute_df[] <- Map(paste, names(attribute_df), attribute_df, sep = '=')
attribute_df<-attribute_df%>%
  tidyr::unite("attributes",ID:wacholder2023_qval, sep = ";", remove = TRUE, na.rm = FALSE)
orf_info%>%as_tibble()%>%
  mutate(type="cds", score=".",
         phase=0)%>%
bind_cols(
  attribute_df
) %>%select(seqid = chromosome, 
                       source,
                       type,
                       start = coor1,
                       end=coor2,
                       score,
                       strand,
                       phase,
                       attributes
                       
)%>%
  arrange(seqid,start)%>%
  
  readr::write_tsv("~/coexpression/data/interim/orf_info.gff")
