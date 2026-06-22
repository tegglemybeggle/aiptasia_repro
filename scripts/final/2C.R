library(DESeq2)
library(topGO)
library(grid)
library(gridExtra)


dds <- readRDS("data/dds.rds")


gaf <- read.delim("genome/go_annotation.gaf", header = FALSE, sep = "\t", quote = "", comment.char = "!", stringsAsFactors = FALSE)


gene_id <- gaf[[3]]
go_id <- gaf[[5]]

geneMap <- tapply(go_id, gene_id, function(x) paste(unique(x), collapse = ","))

writeLines(
  paste(names(geneMap), geneMap, sep = "\t"),
  "data/gene2go.map"
)

gene2GO = readMappings(file = "data/gene2go.map")

geneUniverse <- rownames(dds)

cluster_1_genes <- scan("data/cluster1.txt", what = "")

geneList1 <- factor(as.integer(geneUniverse %in% cluster_1_genes))
names(geneList1) <- geneUniverse

cluster_2_genes <- scan("data/cluster2.txt", what = "")

geneList2 <- factor(as.integer(geneUniverse %in% cluster_2_genes))
names(geneList2) <- geneUniverse



GOdata1 <- new("topGOdata", ontology = "BP", allGenes = geneList1, annot = annFUN.gene2GO, gene2GO = gene2GO)

result1 <- runTest(GOdata1, algorithm = "weight01", statistic = "fisher")

allRes1 <- GenTable(GOdata1, weightFisher = result1, orderBy = "weightFisher", topNodes = 20)

grob1 <- tableGrob(allRes1)
grob1 <- grid.force(grob1)


pdf("figures/2C_1.pdf", height = 0.35 * nrow(allRes1), width = 10)
grid.newpage()
grid.draw(grob1)
dev.off() 

GOdata2 <- new("topGOdata", ontology = "BP", allGenes = geneList2, annot = annFUN.gene2GO, gene2GO = gene2GO)

result2 <- runTest(GOdata2, algorithm = "weight01", statistic = "fisher")

allRes2 <- GenTable(GOdata2, weightFisher = result2, orderBy = "weightFisher", topNodes = 20)

grob2 <- tableGrob(allRes2)
grob2 <- grid.force(grob2)


pdf("figures/2C_2.pdf", height = 0.35 * nrow(allRes2), width = 10)
grid.newpage()
grid.draw(grob2)
dev.off() 


