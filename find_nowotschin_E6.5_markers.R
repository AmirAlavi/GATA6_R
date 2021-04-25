library(Seurat)
library(dplyr)
source("utils.R")

e6.5 <- readRDS("data/Nowotschin_et_al/nowotschin_E6.5_seurat_object.RDS")

# Annotate the cell list based on percentage of mitochondrial genes expressed. 
e6.5[["percent.mt"]] <- PercentageFeatureSet(object = e6.5, pattern = "^mt-")
VlnPlot(object = e6.5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = e6.5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = e6.5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# The following subset is used to remove low-quality cells from the analysis, based on 
# information from the plots created above. 

e6.5 <- subset(x = e6.5, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 6.5)

e6.5 <- NormalizeData(object = e6.5)
e6.5 <- FindVariableFeatures(object = e6.5)
e6.5 <- ScaleData(object = e6.5, features = rownames(x = e6.5))

# The following lines import a list of cell cycle genes that we want to regress out of the
# analysis. Because we are looking for cell types and not cell state, we don't want cells 
# of the same type being clustered in different places because of differences in their cell
# cycle states. We might not do this if we had a particular cell type that was defined by 
# being highly replicative, or vice versa. 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Get mouse orthologs of these genes
m.s.genes <- convertHumanGeneList(s.genes)
m.g2m.genes <- convertHumanGeneList(g2m.genes)


e6.5 <- CellCycleScoring(object = e6.5, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
e6.5 <- ScaleData(object = e6.5, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(x = e6.5))

e6.5 <- RunPCA(object = e6.5, features = VariableFeatures(object = e6.5))


ElbowPlot(object = e6.5)
e6.5 <- JackStraw(object = e6.5, num.replicate = 100)
e6.5 <- ScoreJackStraw(object = e6.5, dims = 1:20)
JackStrawPlot(object = e6.5, dims = 1:20)

e6.5 <- RunTSNE(object = e6.5, dims.use = 1:15)
Idents(e6.5) <- e6.5@meta.data$CellType
DimPlot(object = e6.5, reduction = "tsne")

e6.5.markers <- FindAllMarkers(object = e6.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
e6.5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
saveRDS(e6.5.markers, file = "data/Nowotschin_et_al/E6.5_markers.RDS")
saveRDS(e6.5, file = "data/Nowotschin_et_al/nowotschin_E6.5_seurat_object_processed.RDS")
