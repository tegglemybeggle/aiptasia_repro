library(DESeq2)
library(ComplexHeatmap)
library(circlize)


rlog <- readRDS("data/rlog.rds")
mat <- assay(rlog)
mat_sym <- mat[, grep("Sym", colnames(mat)), drop = FALSE]
t0_mean_sym <- rowMeans(mat_sym[, 1:3, drop = FALSE])
mat_sym_rel <- mat_sym - t0_mean_sym
cluster1 <- readLines("data/cluster1.txt")
cluster2 <- readLines("data/cluster2.txt")
gene_info <- readRDS("data/gene_info.rds")



col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
cell_height = 4

times = c("0", "3", "12", "24", "48", "96")
labels = c()
for (x in 1:6){
  labels = c(labels, "", times[x], "")
}

legend_param = list(
  title = "L2FC",
  at = c(-4, 0, 4),
  labels = gt_render(c("-3 or lower", "0", "3 or higher")),
  legend_height = unit(5, "cm")
)

bottom_anno <- HeatmapAnnotation(
  time = anno_text(
    labels,
    location = 1,
    just = "top",
    rot = 0,
    height = unit(10, "mm")
  )
)

library(dplyr)

motif_counts <- readRDS("data/motif_counts.rds")

heatmaps = c()
motifs = c("NFkB", "HSF")
ids = c(1, 3)

dir.create(
  "figures/3FG",
  recursive = TRUE,
  showWarnings = FALSE
)

for (j in seq_along(ids)){
  
  df <- motif_counts |>
    filter(motif_alt_id == paste0("MEME-", ids[j]) & n >= 3) |>
    arrange(desc(n))
  
  genes <- df$sequence_name[df$sequence_name %in% rownames(mat)]
  mat_i <- mat_sym_rel[genes, ]

  cluster_col <- case_when(
    genes %in% cluster1 ~ "green",
    genes %in% cluster2 ~ "yellow",
    TRUE ~ "white"
  )

  anno_df <- tibble(gene = genes) %>%
    left_join(
      gene_info %>% select(gene, description),
      by = "gene"
    )


  rows = nrow(mat_i)
  height = cell_height * rows
  
  row_anno <- rowAnnotation(

    cluster = anno_simple(
      cluster_col,
      col = c(
        green = "green",
        yellow = "yellow",
        white = "white"
      ),
      width = unit(6, "mm"),
      gp = gpar(col = "black", lwd = 0.5)
    ),

    annotation = anno_text(
      anno_df$description,
      just = "left",
      location = 0,
      gp = gpar(fontsize = 6)
    ),
    show_annotation_name = FALSE
  )



  heatmap_i <- Heatmap(
    mat_i,
    height = unit(height, "mm"),

    show_row_names = FALSE,
    show_column_names = FALSE,

    right_annotation = row_anno,

    heatmap_legend_param = legend_param,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
  )

  heatmaps <- c(heatmaps, heatmap_i)
  png(
    paste0("figures/3FG/", motifs[j], ".png"),
    width = 3000,
    height = (rows * cell_height + 50) / 25.4 * 300,
    res = 300
  )

  draw(heatmap_i)

  dev.off()
}







