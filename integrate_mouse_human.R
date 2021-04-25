library(Seurat)
library(ggplot2)
library(cowplot)
library(plyr)
library(caret)
source("utils.R")
# alignment

# Notes:
# Human Cluster identities:
#   0, 1, 2 exVE
#   4: Epi
#   5: emVE
#   3: maybe mesoderm, but probably not in E6.5 (more likely to be germs cell like)

# Function to build a seurat object, but allow subsetting of the genes and renaming of the genes
seuratObjectCreationWrapper <- function(counts.RDS = "data/Nowotschin_et_al/sc_endoderm_all_cells_counts.RDS",
                                        meta.RDS = "data/Nowotschin_et_al/sc_endoderm_all_cells_metadata.RDS",
                                        project.name = "NowotschinE6.5",
                                        assay = "RNA",
                                        min.cells = 0,
                                        min.features = 0,
                                        genes.use = NULL,
                                        genes.names = NULL) {
  counts <- readRDS(counts.RDS)
  meta <- readRDS(meta.RDS)
  
  if (!is.null(genes.use)) {
    print(paste("shape before:", paste(dim(counts), collapse = ', ')))
    counts <- counts[genes.use, ]
    print(paste("shape after :", paste(dim(counts), collapse = ', ')))
  }
  if (!is.null(genes.names)) {
    print(paste("gene names before:", paste(row.names(counts)[1:10], collapse = ', ')))
    row.names(counts) <- genes.names
    print(paste("gene names after :", paste(row.names(counts)[1:10], collapse = ', ')))
  }
  CreateSeuratObject(counts = counts, project = project.name, assay = assay, meta.data = meta, min.cells = min.cells, min.features = min.features)
}

# First, load Joshua's processed mmB_D5 data to get:
# - the list of all genes used in this data (HGNC symbols)
# - the cluster identites that Joshua found for the cells
load("data/GATA6_D5/GATA6-mmB-Clustering 7-18.RData")
human_clusters <- Idents(mmB_D5)
GATA6_all_genes <- row.names(mmB_D5)
rm(mmB_D5)
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(g2m.genes)
rm(s.genes)

# Second, load an example of the Nowotschin data to get:
# - the list of all genes used in this data (MGI symbols)
nowotschin <- readRDS("data/Nowotschin_et_al/nowotschin_E6.5_seurat_object_processed.RDS")
nowotschin_all_genes <- row.names(nowotschin)
rm(nowotschin)

# For the cross-speices integration analysis, we need to be using a common set of genes.
# Because we have to different species, we will only look at the orthologs
# We need to get this mapping of orthologous genes between human and mouse (from Ensembl)
HGNC2MGI_map <- getHuman2MouseGeneMapping(GATA6_all_genes)

# (keep only orthologous mouse genes that are actually in the Nowotschin data)
common_genes <- HGNC2MGI_map[is.element(HGNC2MGI_map$MGI.symbol, nowotschin_all_genes), ]

# Note this mapping is not one-to-one, and there can by many-to-one and one-to-many maps
# from HGNC to MGI. We should keep all of these, and just duplicate the genes columns where necessary

# Now, in order to ensure that the Seurat objects for the species have the exact same genes in the
# same order for integration, we need to set up the count matrices that way BEFORE we create the Seurat
# Object (seems to be the safest and recommend route)
human <- Read10X(data.dir = "data/GATA6_D5/mmB_D5_Data/", gene.column = 2)
human <- human[, names(human_clusters)] # Keep only the subset of cells that Joshua has labeled (some low quality cells were filtered prior to labeling process)
human <- human[common_genes$HGNC.symbol, ]
human <- CreateSeuratObject(counts = human, min.features = 200, project = "10X_mmB_D5")
human <- AddMetaData(object = human, metadata="human", col.name = "species")
# Apply the cluster labels that Joshua found to the human data
human <- AddMetaData(object = human,
                     metadata=factor(rep("unlabeled", dim(human)[2]), levels = c(levels(human_clusters), "unlabeled")),
                     col.name = "CellType")
human[["CellType"]][names(human_clusters), "CellType"] <- human_clusters

# Build the mouse Seurat object to have the exact same genes (their orthologs) in the same order
# and use the the HGNC names, but keep the original MGI names in metadata for reference later
mouse <- seuratObjectCreationWrapper(counts.RDS = "data/Nowotschin_et_al/nowotschin_E6.5_counts.RDS",
                                     meta.RDS = "data/Nowotschin_et_al/nowotschin_E6.5_metadata.RDS",
                                     min.features = 200,
                                     genes.use = common_genes$MGI.symbol,
                                     genes.names = common_genes$HGNC.symbol)
mouse <- AddMetaData(object = mouse, metadata="mouse", col.name = "species")
mouse[["RNA"]] <- AddMetaData(object = mouse[["RNA"]], metadata=common_genes$MGI.symbol, col.name = "MGI.symbol")

# Finally, proceed with integration, following https://satijalab.org/seurat/v3.1/integration.html

species.list <- c(human, mouse)

for (i in 1:length(species.list)) {
  species.list[[i]] <- NormalizeData(species.list[[i]], verbose = FALSE)
  species.list[[i]] <- FindVariableFeatures(species.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}


species.anchors <- FindIntegrationAnchors(object.list = species.list, dims = 1:30)

species.integrated <- IntegrateData(anchorset = species.anchors, dims = 1:30)
DefaultAssay(species.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
species.integrated <- ScaleData(species.integrated, verbose = FALSE)
species.integrated <- RunPCA(species.integrated, npcs = 30, verbose = FALSE)
species.integrated <- RunUMAP(species.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(species.integrated, reduction = "umap", group.by = "species")
p2 <- DimPlot(species.integrated, reduction = "umap", group.by = "CellType", label = TRUE, 
              repel = TRUE) #+ NoLegend()
plot_grid(p1, p2)

# Transfer labels from mouse to human:
model_knn <- train(x = species.integrated[["pca"]]@cell.embeddings[species.integrated$species == "mouse", ],
                   y = species.integrated$CellType[species.integrated$species == "mouse"],
                   method = "knn")
preds <- predict(object = model_knn, species.integrated[["pca"]]@cell.embeddings[species.integrated$species == "human", ])
names(preds) <- names(species.integrated$CellType[species.integrated$species == "human"])
species.integrated$ProjectedCellType <- preds

plot.data <- species.integrated@meta.data[species.integrated@meta.data$species == "human" ,]
g <- ggplot(plot.data, aes(x = CellType, fill = ProjectedCellType))
g <- g + geom_bar(position = "fill")
g <- g + labs(x = "Human Cluster", y = "Portion of Classifications", fill = "Projected Mouse Cell Type") + ggtitle("Projected Mouse Cell Types for each Human Cluster", "(via KNN classifer after integration)")
g

saveRDS(species.integrated, file = "data/integrated_human_and_mouseE6.5.RDS")

# species.transfer.anchors <- FindTransferAnchors(reference = species.integrated, query = human, 
#                                                 dims = 1:30)
# predictions <- TransferData(anchorset = species.transfer.anchors, refdata = species.integrated$CellType, 
#                             dims = 1:30)
# human <- AddMetaData(human, metadata = predictions)
