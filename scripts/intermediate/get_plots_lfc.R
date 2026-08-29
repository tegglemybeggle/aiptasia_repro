library(DESeq2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(yaml)

gene_info <- readRDS("data/gene_info.rds")
config <- read_yaml("config/config.yaml")
all_results <- readRDS("data/all_results.rds")

timepoints <- c(0, 3, 12, 24, 48, 96)

lfc_data <- tibble()

for (t in timepoints) {

  lfc_data <-
    bind_rows(
      lfc_data,
      all_results[[paste0("SYM_APO_", t, "h")]] |>
        as.data.frame() |>
        rownames_to_column("gene") |>
        transmute(
          gene,
          time = t,
          foldchange = 0 - log2FoldChange
        ) |>
    left_join(
      select(gene_info, gene, description),
      by = "gene")

    )
}


plot_genes <- function(gene_list, title) {
  ggplot(lfc_data %>% filter(gene %in% gene_list),
         aes(x=time, y=foldchange, color=description)) +
  geom_point() +
  geom_line(linewidth=1) +
  scale_x_continuous(breaks = c(0, 3, 12, 24, 48, 96)) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
  labs(y = "Log2 fold difference (Sym/Apo)", x = "Time since initiation of heat stress (h)", title = title)
}

    
