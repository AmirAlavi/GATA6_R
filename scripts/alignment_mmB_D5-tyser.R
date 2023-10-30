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
  print("preprocessing reference...")
  ref <- NormAndHVG(ref)
  print("preprocessing query...")
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
mmB_D5 <- readRDS(snakemake@input[[1]])
# Uncomment next line if you want to analyze cluster 3 as subclusters 3.1 and 3.2
# Idents(mmB_D5) <- mmB_D5$subclusters

# Human Gastrula data (reference)
# tyser <- readRDS(snakemake@input[[2]])
tyser_expr <- readRDS(snakemake@input[[2]])
tyser_labels <- readRDS(snakemake@input[[3]])
tyser_expr <- as.data.frame(t(as.matrix(tyser_expr)))
colnames(tyser_expr) <- tyser_labels$cell_name
rownames(tyser_labels) <- tyser_labels$cell_name
# Create a seurat object 
tyser <- CreateSeuratObject(counts = tyser_expr, project = "Human Gastrula", meta.data = tyser_labels)
Idents(tyser) <- tyser$sub_cluster

refdata <- list(as.character(tyser$sub_cluster))
# refdata <- list(as.character(hu_gas$cluster_id))
names(refdata) <- list("Cell_type")


# cell type/cluster palette
#shifter <- function(x, n = 1) {
#  if (n == 0) x else c(tail(x, -n), head(x, n))
#}
#pal <- DiscretePalette(26, "alphabet2")
#pal <- shifter(pal, 9)
#all.ids <- c(levels(Idents(mmB_D5)), levels(tyser$sub_cluster))
#n.unused.colors <- length(pal) - length(all.ids)
#all.ids <- c(all.ids, rep("none", n.unused.colors))
#names(pal) <- all.ids

all.ids <- c(levels(Idents(mmB_D5)), levels(Idents(tyser)))
pal <- DiscretePalette(length(all.ids))
names(pal) <- all.ids

# *** Do the transfer WITH integration: *******************
combined <- TransferAnnotations(tyser, mmB_D5,
                                refdata,
                                c("tyser", "gata"), "AllCells", TRUE)

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
p2 <- DimPlot(combined, reduction = reduction, label = FALSE, cols = pal) + ggtitle("Original Cell type")# + NoLegend()
# Color by transfered cell types
p3 <- DimPlot(combined, reduction = reduction, group.by = "Cell_type", label = FALSE, cols = pal) + ggtitle("Transfered Cell type")# + NoLegend()

patchwork <- (p1 + p2 + p3)
patchwork <- patchwork + plot_annotation(
  title = "Transfer Tyser et al. Annotations onto GATA6 Data",
  subtitle = "Using Seurat, with integration",
  caption = "We used the Tyser Human Gastrula scRNA-seq count data (1195 cells) as the labeled reference dataset, and used our GATA6 cells (4069) as a query dataset.\nThen we used Seurat's FindTransferAnchors and TransferData functions to annotate the query dataset with Cell Type and Age labels from the reference.",
  tag_levels = 'A'
) & theme(plot.tag = element_text(size = 10))
ggsave(snakemake@output[[1]], width = 20, height = 6, units = "in")

# *** Create classification proportion bar plot
plot.data <- combined@meta.data[combined$orig.ident == "10X_mmB_D5", ]
# g <- ggplot(plot.data, aes(x = seurat_clusters, fill = Cell_type))
g <- ggplot(plot.data, aes(x = subclusters, fill = Cell_type))
g <- g + geom_bar(position = "fill")
g <- g + scale_fill_manual(values = pal)
g <- g + labs(x = "Human Cluster", y = "Portion of Classifications",
              fill = "Projected Tyser et al. Cell Type") +
  ggtitle("Projected Tyser et al. Cell Types for each Human Cluster",
          "(via Seurat Transfer after Integration)")

# ggsave("human_gastrula_subclusters_integrated_transfer_label_distribution.png", width = 8, height = 6, units="in")
ggsave(snakemake@output[[2]], width = 8, height = 6, units="in")
# *********************************************************
