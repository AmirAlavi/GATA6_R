library(Seurat)
library(dplyr)

tyser_expr <- readRDS(snakemake@input[[1]])
tyser_labels <- readRDS(snakemake@input[[2]])
tyser_expr <- as.data.frame(t(as.matrix(tyser_expr)))
colnames(tyser_expr) <- tyser_labels$cell_name
rownames(tyser_labels) <- tyser_labels$cell_name
# Create a seurat object 
tyser <- CreateSeuratObject(counts = tyser_expr, project = "Human Gastrula", meta.data = tyser_labels)
Idents(tyser) <- tyser$sub_cluster

# Annotate the cell list based on percentage of mitochondrial genes expressed. 
tyser[["percent.mt"]] <- PercentageFeatureSet(object = tyser, pattern = "^MT-")
VlnPlot(object = tyser, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = tyser, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = tyser, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# The following subset is used to remove low-quality cells from the analysis, based on 
# information from the plots created above. 

# tyser <- subset(x = tyser, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 6.5)

tyser <- NormalizeData(object = tyser)
tyser <- FindVariableFeatures(object = tyser)
tyser <- ScaleData(object = tyser, features = rownames(x = tyser))

# The following lines import a list of cell cycle genes that we want to regress out of the
# analysis. Because we are looking for cell types and not cell state, we don't want cells 
# of the same type being clustered in different places because of differences in their cell
# cycle states. We might not do this if we had a particular cell type that was defined by 
# being highly replicative, or vice versa. 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tyser <- CellCycleScoring(object = tyser, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tyser <- ScaleData(object = tyser, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(x = tyser))

tyser <- RunPCA(object = tyser, features = VariableFeatures(object = tyser), verbose = FALSE)


ElbowPlot(object = tyser)
tyser <- JackStraw(object = tyser, num.replicate = 100)
tyser <- ScoreJackStraw(object = tyser, dims = 1:20)
JackStrawPlot(object = tyser, dims = 1:20)

tyser <- RunTSNE(object = tyser, dims.use = 1:19)
DimPlot(object = tyser, reduction = "tsne")

tyser.markers <- FindAllMarkers(object = tyser, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tyser.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
tyser.markers <- tyser.markers[tyser.markers$p_val_adj < 0.05, ]
saveRDS(tyser.markers, file = snakemake@output[[1]])
saveRDS(tyser, file = snakemake@output[[2]])
