library(DESeq2)


all_results <- readRDS("data/all_results.rds")

res <- all_results[["SYM_APO_0h"]]


cluster1 <- readLines("data/cluster1.txt")
cluster2 <- readLines("data/cluster2.txt")


library(ggplot2)
library(dplyr)
library(tibble)

df <- data.frame(gene = rownames(res), l2fc = res$log2FoldChange)

df <- df |>
  mutate(l2fc = 0 - l2fc) |>
  mutate(cluster = case_when(
    gene %in% cluster1 ~ "1",
    gene %in% cluster2 ~ "2",
    TRUE                ~ "0"
  ))

p <- ggplot(df, aes(x = l2fc)) +
  geom_density(
  aes(color = "All"),
  linewidth = 1
  ) +
  geom_density(
    data = subset(df, cluster == 1),
    aes(color = "Cluster 1"),
    linewidth = 1
  ) +
  geom_density(
    data = subset(df, cluster == 2),
    aes(color = "Cluster 2"),
    linewidth = 1
    ) +
  scale_color_manual(
    values =  c(
      "All" = "black",
      "Cluster 1" = "yellow",
      "Cluster 2" = "green"
    )
  ) +

  coord_cartesian(xlim = c(-3, 3)) +
  theme(legend.position = "right") +
  labs(x = expression(Log[2] * " Fold-difference Sym/Apo at t = 0"))

pdf("figures/2G.pdf")

print(p)

dev.off()


