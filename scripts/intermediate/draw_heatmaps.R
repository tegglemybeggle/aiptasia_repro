library(ComplexHeatmap)
library(circlize)
library(magick)

dpi = 250
cell_height = 0.5

col_fun <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
dummy_col_fun <- colorRamp2(c(-4, 0, 4), c("white", "white", "white"))

times = c("0", "3", "12", "24", "48", "96")
labels = c()
for (x in 1:6){
  labels = c(labels, "", times[x], "")
}

legend_param = list(
  title = "L2FC",
  at = c(-4, -2, 0, 2, 4),
  labels = gt_render(c("-4 or lower", "-2", "0", "2", "4 or higher")),
  legend_height = unit(5, "cm")
)

dummy_legend = Legend(
  col_fun = dummy_col_fun,
  at = c(-4, -2, 0, 2, 4),
  title = "L2FC",
  labels = gt_render(c("-4 or lower", "-2", "0", "2", "4 or higher")),
  legend_gp = gpar(fill = "white", col = "white"),
  labels_gp = gpar(col = "white"),
  title_gp = gpar(col = "white")
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

draw_heatmaps <- function(s1, a1, s2, a2, filename) {

  n1 = nrow(s1)
  n2 = nrow(s2)

  h1 = cell_height * n1
  h2 = cell_height * n2

  ht_sym1 <- Heatmap(
    s1,
    height = unit(h1, "mm"),

    show_row_names = FALSE,
    show_column_names = FALSE,

    column_title = "Sym",
    row_title = "I",       
    column_title_gp = gpar(fontsize = 28, fontface = "bold"),
    row_title_gp    = gpar(fontsize = 32, fontface = "bold"),
    row_title_rot = 0,

    heatmap_legend_param = legend_param,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
  )

  ht_apo1 <- Heatmap(
    a1,
    height = unit(h1, "mm"),

    show_row_names = FALSE,
    show_column_names = FALSE,

    column_title = "Apo",

    column_title_gp = gpar(fontsize = 28, fontface = "bold"),

    show_heatmap_legend = FALSE,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
  )

  ht_sym2 <- Heatmap(
    s2,
    height = unit(h2, "mm"),

    show_row_names = FALSE,
    show_column_names = FALSE,

    row_title = "II",       
    row_title_gp = gpar(fontsize = 32, fontface = "bold"),
    row_title_rot = 0,

    bottom_annotation = bottom_anno,

    show_heatmap_legend = FALSE,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
  )

  ht_apo2 <- Heatmap(
    a2,
    height = unit(h2, "mm"),

    show_row_names = FALSE,
    show_column_names = FALSE,

    bottom_annotation = bottom_anno,

    show_heatmap_legend = FALSE,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE
  )



  top_row <- ht_sym1 + ht_apo1
  bottom_row <- ht_sym2 + ht_apo2

  pdf(NULL)

  top_drawn <- draw(top_row, heatmap_legend_side = "left")
  bottom_drawn <- draw(bottom_row, heatmap_legend_list = list(dummy_legend), heatmap_legend_side = "left")

  
  h1_in <- grid::convertHeight(ComplexHeatmap:::height(top_drawn), "in", TRUE)
  h2_in <- grid::convertHeight(ComplexHeatmap:::height(bottom_drawn), "in", TRUE)
  w_in <- grid::convertWidth(ComplexHeatmap:::width(bottom_drawn), "in", TRUE)
  
  png("toprow.png", width = w_in * dpi, height = h1_in * dpi, res = dpi )
  draw(top_row, heatmap_legend_side = "left", align_heatmap_legend = "heatmap_top")
  dev.off()

  png("bottomrow.png", width = w_in * dpi, height = h2_in * dpi, res = dpi)
  draw(bottom_row, heatmap_legend_list = list(dummy_legend), heatmap_legend_side = "left", align_heatmap_legend = "heatmap_top")
  dev.off()

  img1 <- image_read("toprow.png")
  img2 <- image_read("bottomrow.png")
  file.remove("toprow.png")
  file.remove("bottomrow.png")

  img_vector <- c(img1, img2)

  final <- image_append(img_vector, stack = TRUE)

  image_write(final, path = filename, format = "png")
}
