library(readr)
library(dplyr)
library(purrr)
library(ggplot2)


motif_counts <- read_tsv("data/fimo_out/fimo.tsv") |>
  mutate(sequence_name = sub("\\(.*", "", sequence_name)) |>
  group_by(sequence_name, motif_alt_id) |>
  summarise(n = n())

saveRDS(motif_counts, "data/motif_counts.rds")
