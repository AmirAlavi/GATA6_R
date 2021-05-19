# Libraries
library(dplyr)
library(Seurat)

xiang <- readRDS(snakemake@input[[1]])

new.annotations <- read.csv(snakemake@input[[2]], row.names = 1)
xiang <- AddMetaData(xiang, metadata = new.annotations)
xiang$For.Analysis <- as.factor(xiang$For.Analysis)
# Remove "Ambiguous" cells from the analysis
Idents(xiang) <- xiang$For.Analysis
xiang <- subset(x = xiang, idents = "Ambiguous", invert = TRUE)

# Annotate the cell list based on percentage of mitochondrial genes expressed. 
xiang[["percent.mt"]] <- PercentageFeatureSet(object = xiang, pattern = "^MT-")
VlnPlot(object = xiang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = xiang, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = xiang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# The following subset is used to remove low-quality cells from the analysis, based on 
# information from the plots created above. 

xiang <- NormalizeData(object = xiang)
xiang <- FindVariableFeatures(object = xiang)
xiang <- ScaleData(object = xiang)

# The following lines import a list of cell cycle genes that we want to regress out of the
# analysis. Because we are looking for cell types and not cell state, we don't want cells 
# of the same type being clustered in different places because of differences in their cell
# cycle states. We might not do this if we had a particular cell type that was defined by 
# being highly replicative, or vice versa. 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# This regression is very expensive to run...
xiang <- CellCycleScoring(object = xiang, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
xiang <- ScaleData(object = xiang, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(x = xiang))

xiang <- RunPCA(object = xiang, features = VariableFeatures(object = xiang))
ElbowPlot(object = xiang)
xiang <- JackStraw(object = xiang, num.replicate = 100)
xiang <- ScoreJackStraw(object = xiang, dims = 1:20)
JackStrawPlot(object = xiang, dims = 1:20)

xiang <- FindNeighbors(xiang, dims = 1:20)
xiang <- RunUMAP(object = xiang, dims = 1:20)

Idents(xiang) <- xiang$For.Analysis
DimPlot(object = xiang, reduction = "umap")

xiang.markers <- FindAllMarkers(object = xiang, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
xiang.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

xiang.markers <- subset(xiang.markers, p_val_adj < 0.05)
saveRDS(xiang.markers, file = snakemake@output[[1]])
saveRDS(xiang, file = snakemake@output[[2]])
