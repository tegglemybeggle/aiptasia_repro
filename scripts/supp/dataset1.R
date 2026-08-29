library(DESeq2)
library(tidyverse)

cluster1 = readLines("data/cluster1.txt")
cluster2 = readLines("data/cluster2.txt")

gene_info = readRDS("data/gene_info.rds")

motif_counts = readRDS("data/motif_counts.rds")

all_results = readRDS("data/all_results.rds")

genes = c(cluster1, cluster2)

motif_wide <- motif_counts |>
  pivot_wider(
    names_from = motif_alt_id,
    values_from = n,
    values_fill = 0
  ) |>
  rename(
    gene = sequence_name,
    NFkB_Sites = `MEME-1`,
    HSF_Sites = `MEME-3`
  )

dataset1 <- tibble(gene = genes)

dataset1 <- dataset1 |>
  left_join(
    select(gene_info, gene, description),
    by = "gene"
  ) |>
  mutate(cluster = gene %in% cluster1 + 1
  ) |>
  left_join(
    select(motif_wide, gene, NFkB_Sites, HSF_Sites),
    by = "gene"
  ) |>
  mutate(
    NFkB_Sites = replace_na(NFkB_Sites, 0)
  ) |>
  mutate(
    HSF_Sites = replace_na(HSF_Sites, 0)
  ) 

for (pop in c("SYM", "APO")){
  for (time in c("3", "12", "24", "48", "96")){
    
    res_name <- paste0(pop, "_0_", time, "h")

    res <- as.data.frame(all_results[[res_name]])

    res <- res |>
      rownames_to_column("gene")
    
    dataset1 <- dataset1 |>
      left_join(
        res |>
          select(
            gene,
            !!paste0("Log2FC_", res_name) := log2FoldChange,
            !!paste0("padj_", res_name) := padj
          ),
        by = "gene"
      )
  }
}

for (time in c("0", "3")){
  
  res_name <- paste0("SYM_APO_", time, "h")
  res <- as.data.frame(all_results[[res_name]])
  res <- res |>
    rownames_to_column("gene")

  dataset1 <- dataset1 |>
      left_join(
        res |>
          select(
            gene,
            !!paste0("Log2FC_", res_name) := log2FoldChange,
            !!paste0("padj_", res_name) := padj
          ),
        by = "gene"
      )
}

dataset1 <- dataset1 |>
  arrange(cluster, desc(Log2FC_SYM_0_3h))


sheet3_genes <- readLines("data/supp/sheet3_genes.txt")

dataset3 <- dataset1 |>
  select(gene, description, cluster, contains("SYM_0_3h"), contains("APO_0_3h"), matches("SYM_APO_(0|3)h$")) |>
  filter(gene %in% sheet3_genes) |>
  arrange(cluster, desc(Log2FC_APO_0_3h))


sheet5_genes <- readLines("data/supp/sheet5_genes.txt")

dataset5 <- dataset1 |>
  select(gene, description, Log2FC_SYM_APO_0h, padj_SYM_APO_0h, cluster) |>
  mutate(`Log2 fold-difference (Sym/Apo 0 h)` = 0 - Log2FC_SYM_APO_0h, .before = Log2FC_SYM_APO_0h) |>
  select(-Log2FC_SYM_APO_0h) |>
  filter(gene %in% sheet5_genes) |>
  filter(`Log2 fold-difference (Sym/Apo 0 h)` >= 0) |>
  arrange(cluster, desc(`Log2 fold-difference (Sym/Apo 0 h)`)) |>
  group_by(cluster) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  relocate(rank, .before = 1)
  
sym_up_genes <- readLines("data/sym_up_genes.txt")
sym_down_genes <- readLines("data/sym_down_genes.txt")
dds <- readRDS("data/dds.rds")

lfc12 <- lfcShrink(dds = dds, coef = "time_12_vs_0", type = "apeglm")

dataset6 <- dataset1 |>
  select(gene, description)

lfc12 <- as.data.frame(lfc12) |> 
  rownames_to_column("gene") |>
  mutate(`Fold-decrease between 0 and 12 h` = 2^(0 - log2FoldChange)) |>
  mutate(`Fold-increase between 0 and 12 h` = 2^(log2FoldChange))

dataset6 <- dataset6 |>
  left_join(lfc12 |>
    select(gene, `Fold-decrease between 0 and 12 h`, padj),
  by = "gene") |>
  filter(gene %in% sym_up_genes & `Fold-decrease between 0 and 12 h` >= 3) |>
  arrange(desc(`Fold-decrease between 0 and 12 h`)) |>
  mutate(rank = row_number(), .before = 1)

dataset7 <- dataset1 |>
  select(gene, description) |>
  left_join(lfc12 |>
    select(gene, `Fold-increase between 0 and 12 h`, padj),
  by = "gene") |>
  filter(gene %in% sym_down_genes & `Fold-increase between 0 and 12 h` >= 3) |>
  arrange(desc(`Fold-increase between 0 and 12 h`)) |>
  mutate(rank = row_number(), .before = 1)


