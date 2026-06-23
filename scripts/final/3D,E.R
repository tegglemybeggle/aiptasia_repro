library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(DESeq2)

all_results <- readRDS("data/all_results.rds")
motif_counts <- readRDS("data/motif_counts.rds")

res <- all_results[["SYM_0_3h"]]

res_df <- data.frame(sequence_name = rownames(res), l2fc = res$log2FoldChange)

motifs <- c("NFkB", "HSF")
ids <- c(1, 3)

plots <- list()

dir.create(
  "figures/3DE",
  recursive = TRUE,
  showWarnings = FALSE
)

for (idx in 1:length(ids)){

  motif_count <- motif_counts |> filter(motif_alt_id == paste0("MEME-", ids[idx]))
 
  df <- left_join(res_df, motif_count, by = "sequence_name")

  df$n[is.na(df$n)] <- 0

  plot_df <- map_dfr(1:5, \(i) {

    df |>
      filter(n >= i) |>
      mutate(label = paste0("N >= ", i))
  })

  p <- ggplot(df, aes(x = l2fc)) + 
    geom_density(aes(color = "All Genes"), linewidth = 1) +
    geom_density(data = plot_df, aes(color = label), linewidth = 1) +
    labs(color = paste0("Putative ", motifs[idx], " Binding Sites"), x = "Log2 Fold Change (3h : 0h)") +
    theme_classic()

  ggsave(
    filename = paste0("figures/3DE/motif_density_", motifs[idx], ".pdf"),
    plot = p,
    width = 5,
    height = 4
  )

}
