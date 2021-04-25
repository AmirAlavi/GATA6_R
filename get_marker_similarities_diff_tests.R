library(Seurat)

source("similarity.R")
source("utils.R")
source("find_human_gastrula_markers_functional.R")


create_marker_list <- function(markers.df) {
  markers <- vector(mode = "list", length = length(unique(markers.df$cluster)))
  names(markers) <- unique(markers.df$cluster)
  for (i in 1:length(markers)) {
    cluster <- names(markers)[i]
    markers[[i]] <- markers.df$gene[markers.df$cluster == cluster]
  }
  return(markers)
}

find_markers <- function(seurat_obj, test = "wilcox", only.pos = TRUE) {
  markers <- FindAllMarkers(object = seurat_obj, test.use = test, only.pos = only.pos, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  markers <- markers[markers$p_val_adj < 0.05, ]
  return(markers)
}

find_similarities <- function(gata_obj, hu_gas_obj, test, only.pos) {
  gata6.markers <- find_markers(gata_obj, test, only.pos)
  gata6.markers <- create_marker_list(gata6.markers)
  gata6.all.genes <- row.names(gata_obj)
  
  hu_gas.markers <- find_markers(hu_gas_obj, test, only.pos)
  hu_gas.markers <- create_marker_list(hu_gas.markers)
  hu_gas.all.genes <- row.names(hu_gas_obj)
  
  background <- intersect(gata6.all.genes, hu_gas.all.genes)
  
  g <- compareMarkersBetweenDatasets(gata6.markers, hu_gas.markers, background, "GATA6_cluster", "Human_Gastrula_cluster")
  g <- g + labs(x = "Human Gastrula cluster", y = "-log10(Hyperg Pval)")
  fname.base <- paste0("2020_12_22/hyperg_human_gastrula_", test)
  ggsave(paste0(fname.base, ".pdf"), plot = g, width = 20, height = 5, units = "in")
  ggsave(paste0(fname.base, ".png"), plot = g, width = 20, height = 5, units = "in")
}

get_datas <- function() {
  # get all possible human genes from our GATA6 experiment:
  load("data/GATA6_D5/GATA6-mmB-Clustering 7-18.RData")
  # Also load the GATA6 markers and put them in a list
  subclusters <- readRDS("data/GATA6_D5/subcluster_3_idents.rds")
  mmB_D5$subcluster_three <- subclusters
  Idents(mmB_D5) <- subclusters
  rm(mmB_D5.data)
  rm(mmB_D5.markers)
  rm(s.genes)
  rm(g2m.genes)
  
  # Human Gastrula
  hu_gas <- load_data()
  hu_gas <- preprocess_data(hu_gas)
  return(list(mmB_D5, hu_gas))
}


datas <- get_datas()
gata6 <- datas[[1]]
hu_gas <- datas[[2]]

tests <- list("wilcox", "bimod", "t", "LR", "MAST", "roc")
# tests <- list("LR", "MAST", "roc")

only.pos = FALSE
for (test.name in tests) {
  print(test.name)
  find_similarities(gata6, hu_gas, test.name, only.pos)
}

