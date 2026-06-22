library(rtracklayer)
library(dplyr)

gtf <- import("genome/genomic.gtf")

gtf_df <- as.data.frame(gtf)

gene_info <- gtf_df %>%
  filter(type == "gene") %>%
  select(
    gene,
    description,
  ) %>%
  distinct()


saveRDS(gene_info, "data/gene_info.rds")


