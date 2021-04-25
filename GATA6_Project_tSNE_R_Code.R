#Seurat tSNE Code for GATA6 Project

# setwd('C:/Users/Josh/SKMO Dropbox/Joshua Hislop/GATA6 Embryogenesis Project/From JH/Single Cell')
# options("BioC_mirror")
# install.packages('Seurat')
library(Seurat)
library(dplyr)
library(Matrix)

# mmB_D5.data = Read10X(data.dir = "Single Cell/Old mmB Data/", gene.column = 2) #Gene column 2 is symbols; setting this number to 1 results in ensembl IDs
mmB_D5.data = Read10X(data.dir = "mmB D5 Data/", gene.column = 2)
mmB_D5 = CreateSeuratObject(counts = mmB_D5.data, min.cells = 5, min.features = 200, project = "10X_mmB_D5")

# Annotate the cell list based on percentage of mitochondrial genes expressed. 

mmB_D5[["percent.mt"]] <- PercentageFeatureSet(object = mmB_D5, pattern = "^MT-")
VlnPlot(object = mmB_D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(object = mmB_D5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = mmB_D5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# The following subset is used to remove low-quality cells from the analysis, based on 
# information from the plots created above. 

mmB_D5 <- subset(x = mmB_D5, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 6.5)

mmB_D5 <- NormalizeData(object = mmB_D5)
mmB_D5 <- FindVariableFeatures(object = mmB_D5)
mmB_D5 <- ScaleData(object = mmB_D5, features = rownames(x = mmB_D5))

# The following lines import a list of cell cycle genes that we want to regress out of the
# analysis. Because we are looking for cell types and not cell state, we don't want cells 
# of the same type being clustered in different places because of differences in their cell
# cycle states. We might not do this if we had a particular cell type that was defined by 
# being highly replicative, or vice versa. 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mmB_D5 <- CellCycleScoring(object = mmB_D5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mmB_D5 <- ScaleData(object = mmB_D5, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(x = mmB_D5))

mmB_D5 <- RunPCA(object = mmB_D5, features = VariableFeatures(object = mmB_D5))

# I normally use the elbow plot to determine the number of dimensions to use in the following 
# "FindNeighbors" function. Because the elbow plot was not extremely clear to me, I followed up 
# with the Jackstraw plot to arrive at 14 dimensions. 

ElbowPlot(object = mmB_D5)
mmB_D5 <- JackStraw(object = mmB_D5, num.replicate = 100)
#mmB_D5 <- ScoreJackStraw(object = mmB_D5, dims = 1:30)
mmB_D5 <- ScoreJackStraw(object = mmB_D5, dims = 1:20)
#JackStrawPlot(object = mmB_D5, dims = 1:30)
JackStrawPlot(object = mmB_D5, dims = 1:20)

mmB_D5 <- FindNeighbors(object = mmB_D5, dims = 1:14)
mmB_D5 <- FindClusters(object = mmB_D5, resolution = 0.5)

mmB_D5 <- RunTSNE(object = mmB_D5, dims.use = 1:14)
DimPlot(object = mmB_D5, reduction = "tsne")

mmB_D5.markers <- FindAllMarkers(object = mmB_D5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mmB_D5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
