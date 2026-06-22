library(DESeq2)


dds <- readRDS("data/dds.rds")

norm_counts <- readRDS("data/rlog.rds")
mat <- assay(norm_counts)

genes <- scan("data/early_genes.txt", what = character())
mat <- mat[genes, , drop = FALSE]

time_order <- order(dds$time)
mat <- mat[, time_order]

mat_sym <- mat[, grep("Sym", colnames(mat)), drop = FALSE]

t0_mean_sym <- rowMeans(mat_sym[, 1:3, drop = FALSE])

mat_sym_rel <- mat_sym - t0_mean_sym

km <- kmeans(mat_sym_rel, 2, nstart=50)


cluster1_genes = rownames(mat_sym_rel)[km$cluster == 1]
cluster2_genes = rownames(mat_sym_rel)[km$cluster == 2]

writeLines(cluster1_genes, "data/cluster1.txt")
writeLines(cluster2_genes, "data/cluster2.txt")


