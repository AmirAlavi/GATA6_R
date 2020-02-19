library(data.table)
library(Seurat)

setwd("joshua_workspace_data")

load("GATA6-mmB-Clustering 7-18.RData")

DimPlot(object = mmB_D5, reduction = "tsne")

data_to_write_out <- as.data.frame(as.matrix(mmB_D5@assays$RNA@scale.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "scale.data.csv")

data_to_write_out <- as.data.frame(as.matrix(mmB_D5@reductions$pca@cell.embeddings))
fwrite(x = data_to_write_out, row.names = TRUE, file = "pca.csv")

data_to_write_out <- as.data.frame(as.matrix(mmB_D5@reductions$tsne@cell.embeddings))
fwrite(x = data_to_write_out, row.names = TRUE, file = "tsne.csv")

data_to_write_out <- as.data.frame(as.matrix(mmB_D5@active.ident))
fwrite(x = data_to_write_out, row.names = TRUE, file = "active.ident.csv")