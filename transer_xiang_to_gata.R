library(Seurat)

library(ggplot2)
library(cowplot)
library(patchwork)

NormAndHVG <- function(data) {
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  return(data)
}

TransferAnnotations <-
  function(ref, query, refdata.list, add.cell.ids = NULL, project = "SeuratProject") {
    ref <- NormAndHVG(ref)
    query <- NormAndHVG(query)
    transfer.anchors <- FindTransferAnchors(reference = ref, query = query, dims = 1:30)
    for (i in 1:length(refdata.list)) {
      refdata <- refdata.list[[i]]
      name <- names(refdata.list[i])
      annotation.predicted <- TransferData(anchorset = transfer.anchors, refdata = refdata, dims = 1:30)
      
      query <- AddMetaData(query, metadata = annotation.predicted$predicted.id, col.name = name)
    }
    combined <- merge(ref, y = query, add.cell.ids = add.cell.ids, project = project, merge.data = TRUE)
    return(combined)
  }

# Our data (query)
load("data/GATA6_D5/GATA6-mmB-Clustering 7-18.RData")
# human_clusters <- Idents(mmB_D5)
# GATA6_all_genes <- row.names(mmB_D5)
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(g2m.genes)
rm(s.genes)
gc()

# Xiang data (reference)
xiang <- readRDS("data/Xiang_et_al_counts/Xiang_counts_seurat_object.RDS")
Idents(xiang) <- xiang$Cell_type

# # Remove unlabled and low-quality cells from the spleen reference
# Idents(spleen) <- spleen$annotation
# spleen <- subset(x = spleen, idents = unique(spleen$annotation[spleen$annotation != "" & spleen$annotation != "low-quality"]))
# lymph <- GetAllDataForTissue("lymph", lymph_datasets)

refdata <- list(as.character(xiang$Cell_type), as.character(xiang$Age))
names(refdata) <- list("Cell_type", "Age")

combined <- TransferAnnotations(xiang, mmB_D5,
                                refdata,
                                c("xng", "gata"), "AllCells")

combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# Color by dataset
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident") + ggtitle("Dataset ID")
# Color by Age
p2 <- DimPlot(combined, reduction = "umap", group.by = "Age", label = TRUE, repel = TRUE) + ggtitle("Transfered Age")
# Color by original cell types
p3 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("Original Cell type")# + NoLegend()
# Color by transfered cell types
p4 <- DimPlot(combined, reduction = "umap", group.by = "Cell_type", label = TRUE, repel = TRUE) + ggtitle("Transfered Cell type")# + NoLegend()

patchwork <- (p1 + p2) / (p3 + p4)
patchwork + plot_annotation(
  title = "Transfer Xiang Annotations onto GATA6 Data",
  subtitle = "Using Seurat, without dataset (batch) alignment",
  caption = "We used the Xiang et al. scRNA-seq count data (555 cells) as the labeled reference dataset, and used our GATA6 cells 4069 as a query dataset.\nThen we used Seurat's FindTransferAnchors and TransferData functions to annotate the query dataset with Cell Type and Age labels from the reference.",
  tag_levels = 'A'
) & theme(plot.tag = element_text(size = 10))

# Best image dimensions are 3000px x 800px (for pdf use 20in x 5.333in)
