library(DESeq2)


dds <- readRDS("data/dds.rds")

norm_counts <- readRDS("data/rlog.rds")
mat <- assay(norm_counts)


mat_sym <- mat[, grep("Sym", colnames(mat)), drop = FALSE]
mat_apo <- mat[, grep("Apo", colnames(mat)), drop = FALSE]


t0_mean_sym <- rowMeans(mat_sym[, 1:3, drop = FALSE])

mat_sym_rel <- mat_sym - t0_mean_sym
mat_apo_rel <- mat_apo - t0_mean_sym


cluster1_genes <- readLines("data/cluster1.txt")
cluster2_genes <- readLines("data/cluster2.txt")



sym1 <- mat_sym_rel[cluster1_genes, ]
sym2 <- mat_sym_rel[cluster2_genes, ]
apo1 <- mat_apo_rel[cluster1_genes, ]
apo2 <- mat_apo_rel[cluster2_genes, ]

o1 <- order(rowMeans(sym1[, 4 : 6, drop = FALSE]), decreasing = TRUE)
o2 <- order(rowMeans(sym2[, 4 : 6, drop = FALSE]), decreasing = TRUE)

sym1 <- sym1[o1, ]
sym2 <- sym2[o2, ]
apo1 <- apo1[o1, ]
apo2 <- apo2[o2, ]


source("scripts/intermediate/draw_heatmaps.R")

draw_heatmaps(sym1, apo1, sym2, apo2, "figures/2A.png")
