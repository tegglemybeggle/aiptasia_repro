
source("scripts/intermediate/get_plots.R")
gene_info <- readRDS("data/gene_info.rds")

cluster1_samples <-gene_info %>%
  filter(
    gene %in% readLines("data/cluster1.txt"),
    !grepl("^uncharacterized", description, ignore.case = TRUE)
  ) %>%
  slice_sample(n = 6) %>%
  pull(gene)

cluster2_samples <-gene_info %>%
  filter(
    gene %in% readLines("data/cluster2.txt"),
    !grepl("^uncharacterized", description, ignore.case = TRUE)
  ) %>%
  slice_sample(n = 6) %>%
  pull(gene)


pdf("figures/2B_1.pdf")

plot_genes(cluster1_samples, 3, 2)

dev.off()

pdf("figures/2B_2.pdf")

plot_genes(cluster2_samples, 3, 2)

dev.off()
