library(Seurat)
library(dplyr)

source("scripts/similarity.R")
source("scripts/utils.R")


loadMaMarkers <- function() {
  markers <- read.csv(snakemake@input[[3]], colClasses = "character", row.names = 1)
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    markers_list[[i]] <- markers$gene[markers$cluster == cur_cluster]
  }
  return(markers_list)
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

# Ma et al. (NHP data)
Ma_markers <- loadMaMarkers()
Ma_counts <- read.csv(snakemake@input[[2]], row.names = 1)
backgroundMa <- intersect(GATA6_all_genes, row.names(Ma_counts))

g <- compareMarkersBetweenDatasets(GATA6_markers, Ma_markers, backgroundMa, "GATA6_cluster", "Ma_cluster")
g + labs(x = "Ma et al. Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[1]], width = 10, height = 5, units = "in")