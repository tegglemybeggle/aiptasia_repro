library(DESeq2)

dds <- readRDS("data/dds.rds")

res_T3vT0 <- results(dds, name = "time_3_vs_0") 

early_genes <- rownames(res_T3vT0)[
  !is.na(res_T3vT0$padj) &
  res_T3vT0$padj < 0.05 &
  res_T3vT0$log2FoldChange > 2
]

writeLines(early_genes, "data/early_genes.txt")

