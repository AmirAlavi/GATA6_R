library(Seurat)
library(data.table)


# This is a large file, so using data.table instead of data.frame
counts <- fread(file = "data/Xiang_et_al_counts/counts.csv")
counts <- as.matrix(counts, rownames=1)
meta <- read.csv("data/Xiang_et_al_counts/meta.csv", header = TRUE, row.names = 1)

xiang <- CreateSeuratObject(counts = counts, project = "Xiang", assay = "RNA", meta.data = meta)

saveRDS(xiang, file = "data/Xiang_et_al_counts/Xiang_counts_seurat_object.RDS")
