library(DESeq2)
library(dplyr)
library(ggplot2)

dds <- readRDS("data/dds.rds")
all_results <- readRDS("data/all_results.rds")

time_levels <- levels(dds$time)


up = readLines("data/sym_up_genes.txt")
down = readLines("data/sym_down_genes.txt")
genes = c(up, down)



df_list <- lapply(time_levels, function(t) {
  res <- all_results[[paste0("SYM_APO_", t, "h")]]
  
  data.frame(
    Gene = genes,
    time = t,
    log2FC = res[genes, "log2FoldChange", drop = TRUE]
  )
})

df <- do.call(rbind, df_list)

df$time <- factor(df$time, levels = c("0", "3", "12", "24", "48", "96"))

df <- df |>
  filter(log2FC > -10 & log2FC < 10) |>
  mutate(`Upregulated in` = if_else(Gene %in% up, "Sym", "Apo")) |>
  mutate(log2FC = 0 - log2FC)

p <- df |>
  ggplot(aes(x = time, y = log2FC, color = `Upregulated in`)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  annotate(
    "segment",
    x = 1:6 - 0.4,
    xend = 1:6 + 0.4,
    y = -Inf,
    yend = -Inf
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Time since initiation of heat stress (h)",
    y = "Log2 fold-difference between symbiotic and aposymbiotic anemones"
  )

pdf("figures/4A.pdf")

print(p)

dev.off()
