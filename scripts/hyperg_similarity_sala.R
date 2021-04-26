library(Seurat)
library(dplyr)

source("scripts/similarity.R")
source("scripts/utils.R")


markerDF2ListHelper <- function(df) {
  markers_list <- as.list(df)
  for (i in 1:length(markers_list)) {
    x <- markers_list[[i]]
    markers_list[[i]] <- convertMouseGeneList(x[x != ""])
  }
  return(markers_list)
}

loadSalaMarkers <- function() {
  markers <- read.csv("data/Sala_et_al/Sala_markers.csv", colClasses = "character")
  markerDF2ListHelper(markers)
}

# get all possible human genes from our GATA6 experiment:
mmB_D5 <- readRDS(snakemake@input[[1]])
GATA6_all_genes <- row.names(mmB_D5)
# Also load the GATA6 markers and put them in a list
Idents(mmB_D5) <- mmB_D5$subclusters
mmB_D5.markers <- FindAllMarkers(object = mmB_D5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mmB_D5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
mmB_D5.markers <- mmB_D5.markers[mmB_D5.markers$p_val_adj < 0.05, ]
n_GATA6_cluster <- length(unique(mmB_D5.markers$cluster))

GATA6_markers <- vector(mode = "list", length = n_GATA6_cluster)
names(GATA6_markers) <- unique(mmB_D5.markers$cluster)
for (i in 1:length(GATA6_markers)) {
  cluster <- names(GATA6_markers)[i]
  GATA6_markers[[i]] <- mmB_D5.markers$gene[mmB_D5.markers$cluster == cluster]
}

# Sala et al.
sala_all_genes <- read.csv(snakemake@input[[2]], header = FALSE, sep = "\t", colClasses = "character", col.names = c("Ensembl", "Symbol"))$Symbol
# create background set in terms of human genes (HGNC)
sala_human_orthologs <- convertMouseGeneList(sala_all_genes)
backgroundSala <- intersect(GATA6_all_genes, sala_human_orthologs)
sala_markers <- loadSalaMarkers()
g <- compareMarkersBetweenDatasets(GATA6_markers, sala_markers, backgroundSala, "GATA6_cluster", "Sala_celltype")
g + labs(x = "Sala et al. Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[1]])