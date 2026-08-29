library(GO.db)
library(AnnotationDbi)
library(yaml)
library(purrr)
library(topGO)
library(tibble)

up_genes <- readLines("data/sym_up_genes.txt")
go_config <- read_yaml("config/go_categories.yaml")
gene2GO <- readMappings("data/gene2go.map")

source("scripts/intermediate/get_plots_lfc.R")

get_descendants <- function(go_id) {

  ontology <- Ontology(go_id)

  descendants <- switch(
    ontology,
    BP = GOBPOFFSPRING[[go_id]],
    MF = GOMFOFFSPRING[[go_id]],
    CC = GOCCOFFSPRING[[go_id]],
    NULL
  )

  descendants <- descendants[!is.na(descendants)]

  unique(c(go_id, descendants))
}

go_config$categories <- map(
  go_config$categories,
  function(category) {
    category$terms <- category$roots |>
      map(get_descendants) |>
      unlist(use.names = FALSE) |>
      unique() |>
      na.omit() |>
      as.character()

    category
  }
)

go_config$categories <- map(
  go_config$categories,
  function(category) {

    category$genes <- c()

    for (gene in up_genes) {

      if (!(gene %in% names(gene2GO))) {
        next
      }

      if (any(gene2GO[[gene]] %in% category$terms)) {
        category$genes <- c(category$genes, gene)
      }
    }

    category
  }
)


set.seed(7)

for (category_name in names(go_config$categories)) {

  category <- go_config$categories[[category_name]]

  go_config$categories[[category_name]]$selected_genes <-
    sample(category$genes, min(4, length(category$genes)))
}


dir.create("figures/4BG")
graphics.off()

for (category in go_config$categories) {

  if (length(category$selected_genes) == 0) {
    next
  }

  pdf(paste0("figures/4BG/", category$label, ".pdf"))
  print(plot_genes(category$selected_genes, category$label))
  dev.off()
}
