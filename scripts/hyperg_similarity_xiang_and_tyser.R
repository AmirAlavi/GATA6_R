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

loadTyserMarkers <- function() {
  markers <- readRDS(snakemake@input[[2]])
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markers$gene[markers$cluster == cur_cluster]
    markers_list[[i]] <- cur_markers
  }
  return(markers_list)
}

loadXiangMarkers <- function() {
  markers <- read.csv("data/Xiang_et_al_counts/Xiang_markers.csv", colClasses = "character")
  markers_list <- as.list(markers)
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


# Tyser et al.
tyser <- readRDS(snakemake@input[[3]])
tyser_all_genes <- row.names(tyser)
background <- intersect(GATA6_all_genes, tyser_all_genes)

tyser_markers <- loadTyserMarkers()

g <- compareMarkersBetweenDatasets(GATA6_markers, tyser_markers, background, "GATA6_cluster", "Tyser_cluster")
g + labs(x = "Tyser et al. cluster", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[1]])

# Xiang et al.
xiang <- readRDS(snakemake@input[[4]])
xiang_all_genes <- row.names(xiang)
background <- intersect(GATA6_all_genes, xiang_all_genes)

xiang_markers <- loadXiangMarkers()

g <- compareMarkersBetweenDatasets(GATA6_markers, xiang_markers, background, "GATA6_cluster", "Xiang_cluster")
g + labs(x = "Xiang et al. cluster", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[2]])

# Tyser et al. and Xiang et al. together
background <- Reduce(intersect, list(GATA6_all_genes, xiang_all_genes, tyser_all_genes))
both_markers <- c(xiang_markers, tyser_markers)
g <- compareMarkersBetweenDatasets(GATA6_markers, both_markers, background, "GATA6_cluster", "Xiang_Tyser_cluster")
g + labs(x = "Xiang et al./Tyser et al. cluster", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[3]])