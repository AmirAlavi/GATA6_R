library(Seurat)

library(ggplot2)
library(cowplot)
library(patchwork)
library(Matrix)


NormAndHVG <- function(data) {
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  return(data)
}

IntegrateDataHelper <- function(object.list) {
  anchors <- FindIntegrationAnchors(object.list = object.list)
  integrated <- IntegrateData(anchorset = anchors)
  return(integrated)
}

TransferAnnotations <- function(ref, query, refdata.list, add.cell.ids = NULL, project = "SeuratProject", integrate = FALSE) {
  ref <- NormAndHVG(ref)
  query <- NormAndHVG(query)
  
  if (integrate) {
    print("integrating data first...")
    integrated <- IntegrateDataHelper(list(ref, query))
    DefaultAssay(integrated) <- "integrated"
    print("subsetting")
    ref <- subset(integrated, cells = colnames(integrated)[integrated$orig.ident == ref$orig.ident[[1]]])
    query <- subset(integrated, cells = colnames(integrated)[integrated$orig.ident == query$orig.ident[[1]]])
    print(dim(query))
    print(dim(ref))
  }
  
  print("transferring...")
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
subclusters <- readRDS("data/GATA6_D5/subcluster_3_idents.rds")
mmB_D5$subcluster_three <- subclusters
Idents(mmB_D5) <- subclusters
# human_clusters <- Idents(mmB_D5)
# GATA6_all_genes <- row.names(mmB_D5)
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(g2m.genes)
rm(s.genes)
gc()

# Human Gastrula data (reference1)
hu_gas_labels <- readRDS("data/Human_Gastrula/umap.rds")
hu_gas_expr <- readRDS("data/Human_Gastrula/express_vals.rds")
hu_gas_expr <- as.data.frame(t(as.matrix(hu_gas_expr)))
colnames(hu_gas_expr) <- hu_gas_labels$cell_name
rownames(hu_gas_labels) <- hu_gas_labels$cell_name
# Create a seurat object 
hu_gas <- CreateSeuratObject(counts = hu_gas_expr, project = "Human Gastrula", meta.data = hu_gas_labels)
Idents(hu_gas) <- hu_gas$cluster_id

# Xiang data (reference2)
xiang <- readRDS("data/Xiang_et_al_counts/Xiang_counts_seurat_object.RDS")
Idents(xiang) <- xiang$Cell_type

# Remove the irrelevant cell-fates
xiang <- subset(xiang, idents = c("EVTs", "CTBs", "STBs"), invert = TRUE)

# # Remove unlabled and low-quality cells from the spleen reference
# Idents(spleen) <- spleen$annotation
# spleen <- subset(x = spleen, idents = unique(spleen$annotation[spleen$annotation != "" & spleen$annotation != "low-quality"]))
# lymph <- GetAllDataForTissue("lymph", lymph_datasets)








# refdata <- list(as.character(hu_gas$sub_cluster))
refdata <- list(as.character(hu_gas$cluster_id))
names(refdata) <- list("Cell_type")

# *** Do the transfer WITHOUT integration: ****************
combined <- TransferAnnotations(hu_gas, mmB_D5,
                                refdata,
                                c("hu_gas", "gata"), "AllCells")

combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# reduction = "umap"
reduction = "pca"
# Color by dataset
p1 <- DimPlot(combined, reduction = reduction, group.by = "orig.ident") + ggtitle("Dataset ID")
# Color by original cell types
p2 <- DimPlot(combined, reduction = reduction, label = TRUE, repel = TRUE) + ggtitle("Original Cell type")# + NoLegend()
# Color by transfered cell types
p3 <- DimPlot(combined, reduction = reduction, group.by = "Cell_type", label = TRUE, repel = TRUE) + ggtitle("Transfered Cell type")# + NoLegend()

patchwork <- (p1 + p2 + p3) 
patchwork + plot_annotation(
  title = "Transfer Human Gastrula Annotations onto GATA6 Data",
  subtitle = "Using Seurat, without dataset (batch) alignment",
  caption = "We used the Human Gastrula scRNA-seq count data (1195 cells) as the labeled reference dataset, and used our GATA6 cells 4069 as a query dataset.\nThen we used Seurat's FindTransferAnchors and TransferData functions to annotate the query dataset with Cell Type labels from the reference.",
  tag_levels = 'A'
) & theme(plot.tag = element_text(size = 10))
ggsave(paste0("human_gastrula_onto_GATA6_", reduction, "_align_off.png"),
       width = 10, height = 4, units="in")
# *********************************************************

# *** Do the transfer WITH integration: *******************
combined <- TransferAnnotations(hu_gas, mmB_D5,
                                refdata,
                                c("hu_gas", "gata"), "AllCells", TRUE)

blank.mat.sparse <- as(matrix(, nrow=0, ncol = 0), "dgCMatrix")
combined[["integrated"]]@counts <- blank.mat.sparse
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE, features = rownames(combined))
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# Color by dataset
reduction = "umap"
# reduction = "pca"
p1 <- DimPlot(combined, reduction = reduction, group.by = "orig.ident") + ggtitle("Dataset ID")
# Color by original cell types
p2 <- DimPlot(combined, reduction = reduction, label = TRUE, repel = TRUE) + ggtitle("Original Cell type")# + NoLegend()
# Color by transfered cell types
p3 <- DimPlot(combined, reduction = reduction, group.by = "Cell_type", label = TRUE, repel = TRUE) + ggtitle("Transfered Cell type")# + NoLegend()

patchwork <- (p1 + p2 + p3)
patchwork + plot_annotation(
  title = "Transfer Human Gastrula Annotations onto GATA6 Data",
  subtitle = "Using Seurat, with integration",
  caption = "We used the Human Gastrula scRNA-seq count data (1195 cells) as the labeled reference dataset, and used our GATA6 cells 4069 as a query dataset.\nThen we used Seurat's FindTransferAnchors and TransferData functions to annotate the query dataset with Cell Type and Age labels from the reference.",
  tag_levels = 'A'
) & theme(plot.tag = element_text(size = 10))

# *** Create classification proportion bar plot
plot.data <- combined@meta.data[combined$orig.ident == "10X_mmB_D5", ]
# g <- ggplot(plot.data, aes(x = seurat_clusters, fill = Cell_type))
g <- ggplot(plot.data, aes(x = subcluster_three, fill = Cell_type))
g <- g + geom_bar(position = "fill")
g <- g + labs(x = "Human Cluster", y = "Portion of Classifications",
              fill = "Projected Human Gastrula Cell Type") +
  ggtitle("Projected Human Gastrula Cell Types for each Human Cluster",
          "(via Seurat Transfer after Integration)")
g
ggsave("human_gastrula_subclusters_integrated_transfer_label_distribution.png", width = 8, height = 6, units="in")
# *********************************************************










refdata <- list(as.character(xiang$Cell_type), as.character(xiang$Age))
names(refdata) <- list("Cell_type", "Age")

# *** Do the transfer WITHOUT integration: ****************
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
# *********************************************************

# *** Do the transfer WITH integration: *******************
combined <- TransferAnnotations(xiang, mmB_D5,
                                refdata,
                                c("xng", "gata"), "AllCells", TRUE)

blank.mat.sparse <- as(matrix(, nrow=0, ncol = 0), "dgCMatrix")
combined[["integrated"]]@counts <- blank.mat.sparse
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE, features = rownames(combined))
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
# Color by dataset
# reduction = "umap"
reduction = "pca"
p1 <- DimPlot(combined, reduction = reduction, group.by = "orig.ident") + ggtitle("Dataset ID")
# Color by Age
p2 <- DimPlot(combined, reduction = reduction, group.by = "Age", label = TRUE, repel = TRUE) + ggtitle("Transfered Age")
# Color by original cell types
p3 <- DimPlot(combined, reduction = reduction, label = TRUE, repel = TRUE) + ggtitle("Original Cell type")# + NoLegend()
# Color by transfered cell types
p4 <- DimPlot(combined, reduction = reduction, group.by = "Cell_type", label = TRUE, repel = TRUE) + ggtitle("Transfered Cell type")# + NoLegend()

patchwork <- (p1 + p2) / (p3 + p4)
patchwork + plot_annotation(
  title = "Transfer Xiang Annotations onto GATA6 Data",
  subtitle = "Using Seurat, with integration",
  caption = "We used the Xiang et al. scRNA-seq count data (555 cells) as the labeled reference dataset, and used our GATA6 cells 4069 as a query dataset.\nThen we used Seurat's FindTransferAnchors and TransferData functions to annotate the query dataset with Cell Type and Age labels from the reference.",
  tag_levels = 'A'
) & theme(plot.tag = element_text(size = 10))

# *** Create classification proportion bar plot
plot.data <- combined@meta.data[combined$orig.ident == "10X_mmB_D5", ]
g <- ggplot(plot.data, aes(x = seurat_clusters, fill = Cell_type))
g <- g + geom_bar(position = "fill")
g <- g + labs(x = "Human Cluster", y = "Portion of Classifications",
              fill = "Projected Xiang Cell Type") +
  ggtitle("Projected Xiang et al. Cell Types for each Human Cluster",
          "(via Seurat Transfer after Integration)")
g
# *********************************************************