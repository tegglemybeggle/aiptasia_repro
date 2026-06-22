library(DESeq2)


dds <- readRDS("data/dds.rds")

norm_counts <- rlog(dds, blind = FALSE)

saveRDS(norm_counts, "data/rlog.rds")
