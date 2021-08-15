library(Seurat)
library(plyr)
load(snakemake@input[[1]])
DimPlot(object = mmB_D5, reduction = "tsne")

# Focus on cluster 3
cluster3 <- mmB_D5[, Idents(mmB_D5) == 3]
DimPlot(object = cluster3, reduction = "tsne")

# Conduct Hierarchical clustering on the scaled data
dist_mat <- dist(t(cluster3@assays$RNA@scale.data), method = "euclidean")
clusters <- hclust(dist_mat, method = "ward.D2")
# Visualize the dendogram
plot(clusters, labels = FALSE)

# Looks like there might be two distinct sub-clusters, pick that clustering
clusterCut <- cutree(clusters, 2)
Idents(object = cluster3) <- clusterCut
DimPlot(object = cluster3, reduction = "tsne")

# Now apply these cluster IDs back to the whole data to find markers
subcluster_idents <- revalue(Idents(object = cluster3), c("1"=3.1, "2"=3.2))
combined_idents <- Idents(mmB_D5)
levels(combined_idents) <- c(levels(combined_idents), levels(subcluster_idents))
combined_idents[attr(subcluster_idents, "names")] <- subcluster_idents
Idents(mmB_D5) <- combined_idents
mmB_D5$subclusters <- combined_idents
DimPlot(object = mmB_D5, reduction = "tsne")
mmB_D5 <- RunUMAP(mmB_D5, dims = 1:20)
DimPlot(object = mmB_D5, reduction = "umap")
saveRDS(mmB_D5, file = snakemake@output[[1]])

# Find markers for the new subclusters
subcluster3.1.markers <- FindMarkers(mmB_D5, ident.1 = 3.1, only.pos = TRUE, min.pct = 0.25)
subcluster3.1.markers.sig <- subcluster3.1.markers[subcluster3.1.markers$p_val_adj < 0.05, ]
dim(subcluster3.1.markers.sig)
head(subcluster3.1.markers.sig, n = 10)

subcluster3.2.markers <- FindMarkers(mmB_D5, ident.1 = 3.2, only.pos = TRUE, min.pct = 0.25)
subcluster3.2.markers.sig <- subcluster3.2.markers[subcluster3.2.markers$p_val_adj < 0.05, ]
dim(subcluster3.2.markers.sig)
head(subcluster3.2.markers.sig, n = 10)

