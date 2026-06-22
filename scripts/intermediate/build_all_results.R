library(DESeq2)
library(dplyr)

dds <- readRDS("data/dds.rds")

time_levels <- levels(dds$time)
res_list <- list()

for (t in time_levels[-1]){

  #Sym vs 0
  res_list[[paste0("SYM_0_", t, "h")]] <- results(dds, name = paste0("time_", t, "_vs_0"))

    #Apo vs 0
  res_list[[paste0("APO_0_", t, "h")]] <- results(dds, contrast = list(c(paste0("time_", t, "_vs_0"), paste0("popApo.time", t))))

}


res_list[["SYM_APO_0h"]] <- results(dds, name = "pop_Apo_vs_Sym")

for (t in time_levels[-1]){

  #Apo vs Sym
  res_list[[paste0("SYM_APO_", t, "h")]] <- results(dds, contrast = list(c("pop_Apo_vs_Sym", paste0("popApo.time", t))))
}


saveRDS(res_list, "data/all_results.rds")


