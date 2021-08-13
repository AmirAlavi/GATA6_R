library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

loadXiangNewMarkers <- function() {
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

PlotModuleScores <- function(obj, markers) {
  for (cell.type in names(markers)) {
    obj <- AddModuleScore(obj, features = markers[cell.type], name = cell.type, search = TRUE)
    names(obj@meta.data)[names(obj@meta.data) == paste0(cell.type, "1")] <- paste0(cell.type, ".score")
  }
  FeaturePlot(obj, reduction = "tsne", features = paste0(names(markers), ".score"), combine = FALSE)
}

new_xiang_markers <- loadXiangNewMarkers()

mmB_D5 <- readRDS(snakemake@input[[1]])
Idents(mmB_D5) <- mmB_D5$subclusters

# cell type/cluster palette
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
pal <- DiscretePalette(26, "alphabet2")
pal <- shifter(pal, 9)
all.ids <- levels(Idents(mmB_D5))
n.unused.colors <- length(pal) - length(all.ids)
all.ids <- c(all.ids, rep("none", n.unused.colors))
names(pal) <- all.ids

DimPlot(mmB_D5, reduction="tsne", cols = pal) + PlotModuleScores(mmB_D5, new_xiang_markers) + plot_annotation(title = "New Xiang Cell type scores")
ggsave(snakemake@output[[1]], scale = 2)
