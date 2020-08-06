source("similarity.R")
source("utils.R")

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
  markers <- readRDS("data/Nowotschin_et_al/E6.5_markers.RDS")
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

loadNHPMarkers <- function() {
  markers <- read.csv("data/NHP/aax7890-Ma-SM-Table-S6.csv", colClasses = "character", row.names = 1)
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    markers_list[[i]] <- markers$gene[markers$cluster == cur_cluster]
  }
  return(markers_list)
}

loadSalaMarkers <- function() {
  markers <- read.csv("data/Sala_et_al/Sala_markers.csv", colClasses = "character")
  markerDF2ListHelper(markers)
}

loadXiangMarkers <- function() {
  markers <- read.csv("data/Xiang_et_al_counts/Xiang_markers.csv", colClasses = "character")
  markers_list <- as.list(markers)
  return(markers_list)
}


# get all possible human genes from our GATA6 experiment:
load("data/GATA6_D5/GATA6-mmB-Clustering 7-18.RData")
GATA6_all_genes <- row.names(mmB_D5)
# Also load the GATA6 markers and put them in a list
n_GATA6_cluster <- length(unique(mmB_D5.markers$cluster))
GATA6_markers <- vector(mode = "list", length = n_GATA6_cluster)
names(GATA6_markers) <- unique(mmB_D5.markers$cluster)
for (i in 1:length(GATA6_markers)) {
  cluster <- names(GATA6_markers)[i]
  GATA6_markers[[i]] <- mmB_D5.markers$gene[mmB_D5.markers$cluster == cluster]
}
rm(mmB_D5)
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(s.genes)
rm(g2m.genes)

# get all possible mouse genes from the mouse endoderm explorer
nowotschin <- readRDS("data/Nowotschin_et_al/nowotschin_E6.5_seurat_object_processed.RDS")
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

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E5.5_markers, background, "GATA6_cluster", "Nowotschin_E5.5_celltype")
g + labs(x = "Nowotschin E5.5 Cell Type", y = "-log10(Hyperg Pval)")

g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E6.5_markers, background, "GATA6_cluster", "Nowotschin_E6.5_celltype")
g + labs(x = "Nowotschin E6.5 Cell Type", y = "-log10(Hyperg Pval)")

# Try with making the marker lists disjoint
nowotschin_E6.5_markers_disjoint <- makeMarkersDisjoint(nowotschin_E6.5_markers)
g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E6.5_markers_disjoint, background, "GATA6_cluster", "Nowotschin_E6.5_celltype")
g + labs(x = "Nowotschin E6.5 Cell Type", y = "-log10(Hyperg Pval)")


g <- compareMarkersBetweenDatasets(GATA6_markers, nowotschin_E7.5_markers, background, "GATA6_cluster", "Nowotschin_E7.5_celltype")
g + labs(x = "Nowotschin E7.5 Cell Type", y = "-log10(Hyperg Pval)")


compareMarkersWithinDataset(GATA6_markers)
compareMarkersWithinDataset(nowotschin_E6.5_markers)
compareMarkersWithinDataset(nowotschin_E6.5_markers_disjoint)

# Xiang et al.
xiang <- readRDS("data/Xiang_et_al_counts/Xiang_counts_seurat_object.RDS")
xiang_all_genes <-row.names(xiang)
background <- intersect(GATA6_all_genes, xiang_all_genes)

xiang_markers <- loadXiangMarkers()

g <- compareMarkersBetweenDatasets(GATA6_markers, xiang_markers, background, "GATA6_cluster", "Xiang_cluster")
g + labs(x = "Xiang et al. cluster", y = "-log10(Hyperg Pval)")


# Sala et al.
sala_all_genes <- read.csv("data/Sala_et_al/genes.tsv", header = FALSE, sep = "\t", colClasses = "character", col.names = c("Ensembl", "Symbol"))$Symbol
# create background set in terms of human genes (HGNC)
sala_human_orthologs <- convertMouseGeneList(sala_all_genes)
backgroundSala <- intersect(GATA6_all_genes, sala_human_orthologs)
sala_markers <- loadSalaMarkers()
g <- compareMarkersBetweenDatasets(GATA6_markers, sala_markers, backgroundSala, "GATA6_cluster", "Sala_celltype")
g + labs(x = "Sala Cell Type", y = "-log10(Hyperg Pval)")

# NHP data
NHP_markers <- loadNHPMarkers()
NHP_counts <- read.csv('data/NHP/Counts Matrix - GSE130114_MF1453.csv', row.names = 1)
backgroundNHP <- intersect(GATA6_all_genes, row.names(NHP_counts))

g <- compareMarkersBetweenDatasets(GATA6_markers, NHP_markers, backgroundNHP, "GATA6_cluster", "NHP_cluster")
g + labs(x = "NHP Cell Type", y = "-log10(Hyperg Pval)")
