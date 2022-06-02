library(Seurat)
library(dplyr)

source("scripts/similarity.R")
source("scripts/utils.R")

loadTyserMarkers <- function() {
  markers <- readRDS(snakemake@input[[5]])
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

markersDF2List <- function(markersDF) {
  n_clusters <- length(unique(markersDF$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markersDF$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markersDF$gene[markersDF$cluster == cur_cluster]
    markers_list[[i]] <- cur_markers
  }
  return(markers_list)
}

xiang <- readRDS(snakemake@input[[1]])
xiang_markers <- readRDS(snakemake@input[[2]])
xiang_all_genes <- row.names(xiang)
iDiscoid_markers <- readRDS(snakemake@input[[3]])
iDiscoid_gene_sets <- readRDS(snakemake@input[[4]])
tyser <- readRDS(snakemake@input[[6]])
tyser_all_genes <- row.names(tyser)
tyser_markers <- loadTyserMarkers()

for (i in 1:length(iDiscoid_markers)) {
  our_day <- names(iDiscoid_markers)[i]
  our_markers <- markersDF2List(iDiscoid_markers[[our_day]])
  our_gene_set <- iDiscoid_gene_sets[[our_day]]

  background <- intersect(our_gene_set, tyser_all_genes)
  g <- compareMarkersBetweenDatasets(our_markers, tyser_markers, background, "iDiscoid_cluster", "Tyser_cluster")
  g + labs(x = "Tyser et al. cluster", y = "-log10(Hyperg Pval)")
  path <- paste0(dirname(snakemake@output[[1]]), "/iDiscoid_", our_day, "_Tyser.png")
  ggsave(path, width = 10, height = 5, units = "in")

  for (j in 1:length(xiang_markers)) {
    xiang_day <- names(xiang_markers)[j]
    xiang_day_markers <- markersDF2List(xiang_markers[[xiang_day]])
    background <- intersect(our_gene_set, xiang_all_genes)
    g <- compareMarkersBetweenDatasets(our_markers, xiang_day_markers, background, "iDiscoid_cluster", "Xiang_cluster")
    g + labs(x = "(new) Xiang et al. cluster", y = "-log10(Hyperg Pval)")
    path <- paste0(dirname(snakemake@output[[1]]), "/iDiscoid_", our_day, "_Xiang_", xiang_day, ".png")
    print(path)
    ggsave(path, width = 10, height = 5, units = "in")
  }
}

file.create(snakemake@output[[1]])