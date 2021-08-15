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

loadNowotschinE4.5Markers <- function() {
  markers <- read.csv("data/Nowotschin_et_al/E4.5_markers.csv", colClasses = "character")
  markerDF2ListHelper(markers)
}

loadNowotschinE5.5Markers <- function() {
  markers <- read.csv("data/Nowotschin_et_al/E5.5_markers.csv", colClasses = "character")
  markerDF2ListHelper(markers)
}

loadNowotschinE6.5Markers <- function() {
  markers <- readRDS(snakemake@input[[3]])
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markers$gene[markers$cluster == cur_cluster]
    markers_list[[i]] <- convertMouseGeneList(cur_markers)
  }
  return(markers_list)
}

loadNowotschinE7.5Markers <- function() {
  markers <- read.csv("data/Nowotschin_et_al/E7.5_markers.csv", colClasses = "character")
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

# get all possible mouse genes from the mouse endoderm explorer
nowotschin <- readRDS(snakemake@input[[2]])
nowotschin_all_genes <- row.names(nowotschin)
rm(nowotschin)

# create background set in terms of human genes (HGNC)
nowotschin_human_orthologs <- convertMouseGeneList(nowotschin_all_genes)
background <- intersect(GATA6_all_genes, nowotschin_human_orthologs)

nowotschin_E4.5_markers <- loadNowotschinE4.5Markers()
nowotschin_E5.5_markers <- loadNowotschinE5.5Markers()
nowotschin_E6.5_markers <- loadNowotschinE6.5Markers()
nowotschin_E7.5_markers <- loadNowotschinE7.5Markers()

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E4.5_markers, background, "GATA6_cluster", "Nowotschin_E4.5_celltype")
g + labs(x = "Nowotschin E4.5 Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[1]], width = 10, height = 5, units = "in")

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E5.5_markers, background, "GATA6_cluster", "Nowotschin_E5.5_celltype")
g + labs(x = "Nowotschin E5.5 Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[2]], width = 10, height = 5, units = "in")

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E6.5_markers, background, "GATA6_cluster", "Nowotschin_E6.5_celltype")
g + labs(x = "Nowotschin E6.5 Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[3]], width = 10, height = 5, units = "in")

# Take in a marker list and make all of the marker sets in it disjoint wrt each other
makeMarkersDisjoint <- function(x) {
  newx <- vector(mode = "list", length = length(x))
  names(newx) <- names(x)
  for (i in 1:length(x)) {
    markers <- x[[i]]
    for (j in 1:length(x)) {
      if (i == j) {
        next
      }
      other_markers <- x[[j]]
      markers <- setdiff(markers, other_markers)
    }
    newx[[i]] <- markers
  }
  return(newx)
}

# Try with making the marker lists disjoint
nowotschin_E6.5_markers_disjoint <- makeMarkersDisjoint(nowotschin_E6.5_markers)
g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E6.5_markers_disjoint, background, "GATA6_cluster", "Nowotschin_E6.5_celltype")
g + labs(x = "Nowotschin E6.5 Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[4]], width = 10, height = 5, units = "in")

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E7.5_markers, background, "GATA6_cluster", "Nowotschin_E7.5_celltype")
g + labs(x = "Nowotschin E7.5 Cell Type", y = "-log10(Hyperg Pval)")
ggsave(snakemake@output[[5]], width = 10, height = 5, units = "in")

compareMarkersWithinDataset(GATA6_markers)
ggsave(snakemake@output[[6]])

compareMarkersWithinDataset(nowotschin_E6.5_markers)
ggsave(snakemake@output[[7]])

compareMarkersWithinDataset(nowotschin_E6.5_markers_disjoint)
ggsave(snakemake@output[[8]])