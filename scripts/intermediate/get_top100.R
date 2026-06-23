library(DESeq2)

dds <- readRDS("data/dds.rds")

lfc <- lfcShrink(dds=dds, coef="time_3_vs_0", type="apeglm")

top100 <- rownames(lfc[order(lfc$log2FoldChange, decreasing=TRUE), ])[1:100]


write.table(top100, file="data/top100_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

