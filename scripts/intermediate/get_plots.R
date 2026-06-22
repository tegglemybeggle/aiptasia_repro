library(DESeq2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

gene_info <- readRDS("data/gene_info.rds")


dds <- readRDS("data/dds.rds")

mat <- counts(dds, normalized=TRUE)

df <- as.data.frame(mat)
df$gene <- rownames(df)

long_df <- df %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  separate(sample, into = c("pop", "time", "rep"), sep = "_") %>%
  mutate(time = as.numeric(time)) %>%
  left_join(
    select(gene_info, gene, description),
    by = "gene")

head(long_df)

long_avgd <- long_df %>%
  group_by(gene, description, pop, time) %>%
  summarize(count = mean(count), .groups= "drop") %>%
  mutate(
    label = paste0(
      str_wrap(description, width = 25),
      "\n",
      gene
    )
  )

head(long_avgd)

plot_genes <- function(gene_list, cols, rows) {
  ggplot(long_avgd %>% filter(gene %in% gene_list),
         aes(x=time, y=count, color=pop, group=pop)) +
  geom_line(linewidth=1) +
  scale_color_manual(values = c("Sym" = "blue", "Apo" = "red")) +
  scale_x_continuous(breaks = c(0, 3, 12, 24, 48, 96)) + 
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  facet_wrap(~label, nrow = rows, ncol = cols, scales = "free_y") 
}

    
