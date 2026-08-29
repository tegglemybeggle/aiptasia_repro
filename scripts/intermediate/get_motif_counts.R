library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

original <- read_tsv("data/fimo_out/fimo.tsv") |>
  filter(motif_alt_id != "MEME_3")

trimmed <- read_tsv("fimo_test_out/fimo.tsv")

fimo <- bind_rows(original, trimmed)

motif_counts <- fimo |>
  mutate(sequence_name = sub("\\(.*", "", sequence_name)) |>
  group_by(sequence_name, motif_alt_id) |>
  summarise(n = n())

saveRDS(motif_counts, "data/motif_counts.rds")
