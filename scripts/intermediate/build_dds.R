library(DESeq2)
library(pheatmap)
library(yaml)

config <- read_yaml("config/config.yaml")
big_data_dir <- config$big_data_dir

directory <- paste0(big_data_dir, "counts")
files <- list.files(directory, pattern = "\\.txt$", full.names = FALSE)
sampleNames <- tools::file_path_sans_ext(files)

ms <- regexec("^(Apo|Sym)_([0-9]+)_[1-3]", sampleNames)
parts <- regmatches(sampleNames, ms)

pop <- vapply(parts, `[[`, character(1), 2)
time <- vapply(parts, `[[`, character(1), 3)

pop <- factor(pop, levels = c("Sym", "Apo"))
time <- factor(time, levels = c("0", "3", "12", "24", "48", "96"))

sampleTable <- data.frame(
  sampleName = sampleNames,
  fileName = files,
  pop = pop,
  time = time,
  stringsAsFactors = FALSE
)

sampleTable <- sampleTable[order(sampleTable$pop, sampleTable$time), ]

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~ pop * time
)

dds <- DESeq(dds)

saveRDS(dds, file = "data/dds.rds")
