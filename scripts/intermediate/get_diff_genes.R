library(DESeq2)

all_results <- readRDS("data/all_results.rds")

res <- all_results[["SYM_APO_0h"]]

sym_up_genes <- rownames(res)[
  !is.na(res$padj) &
  res$padj < 0.05 &
  ((res$log2FoldChange < -2))]

sym_down_genes <- rownames(res)[
  !is.na(res$padj) &
  res$padj < 0.05 &
  ((res$log2FoldChange > 2))]


writeLines(sym_up_genes, "data/sym_up_genes.txt")
writeLines(sym_down_genes, "data/sym_down_genes.txt")
