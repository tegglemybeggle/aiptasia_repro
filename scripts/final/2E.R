library(DESeq2)

data_dir = "~/cleves/new/amid/data/"
dds <- readRDS(paste0(data_dir, "dds.rds"))

norm_counts <- readRDS(paste0(data_dir, "rlog.rds"))
mat <- assay(norm_counts)


mat_sym <- mat[, grep("Sym", colnames(mat)), drop = FALSE]
mat_apo <- mat[, grep("Apo", colnames(mat)), drop = FALSE]


t0_mean_apo <- rowMeans(mat_apo[, 1:3, drop = FALSE])

mat_sym_rel <- mat_sym - t0_mean_apo
mat_apo_rel <- mat_apo - t0_mean_apo


cluster1_genes <- readLines(paste0(data_dir, "cluster1.txt"))
cluster2_genes <- readLines(paste0(data_dir, "cluster2.txt"))


sym1 <- mat_sym_rel[cluster1_genes, ]
sym2 <- mat_sym_rel[cluster2_genes, ]
apo1 <- mat_apo_rel[cluster1_genes, ]
apo2 <- mat_apo_rel[cluster2_genes, ]

cluster1_low <- rowMeans(apo1[, 4:6, drop = FALSE]) <= 1
cluster2_low <- rowMeans(apo2[, 4:6, drop = FALSE]) <= 1

sym1 <- sym1[cluster1_low, ]
sym2 <- sym2[cluster2_low, ]
apo1 <- apo1[cluster1_low, ]
apo2 <- apo2[cluster2_low, ]

o1 <- order(rowMeans(apo1[, 4 : 6, drop = FALSE]), decreasing = TRUE)
o2 <- order(rowMeans(apo2[, 4 : 6, drop = FALSE]), decreasing = TRUE)

sym1 <- sym1[o1, ]
sym2 <- sym2[o2, ]
apo1 <- apo1[o1, ]
apo2 <- apo2[o2, ]


source("~/cleves/new/amid/scripts/get_plots.R")

c1_genes = sample(rownames(sym1), size = 3)
c2_genes = sample(rownames(sym2), size = 3)

pdf("~/cleves/new/figures/2E.pdf")

plot_genes(c(c1_genes, c2_genes), 3, 2) 

dev.off()
