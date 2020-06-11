library(Seurat)
library(data.table)


dt <- fread(file = "data/Xiang_et_al/GSE136447_555-samples-fpkm.txt")
meta <- read.csv("data/Xiang_et_al/41586_2019_1875_MOESM10_ESM.csv", skip = 2, row.names = 1)[, 0:3]
fpkm <- dt[, row.names(meta), with=FALSE]
fpkm <- as.matrix(fpkm, rownames.value=dt$`Gene_ID`)

xiang <- CreateSeuratObject(counts = fpkm, project = "Xiang", assay = "RNA", meta.data = meta)
colnames(xiang[["RNA"]]) <- dt$Gene_ID
xiang[["RNA"]] <- AddMetaData(object = xiang[["RNA"]], metadata = dt$`Gene ID`, col.name = "ensemble.id")
xiang[["RNA"]] <- AddMetaData(object = xiang[["RNA"]], metadata = dt$`Gene Name`, col.name = "symbol")

saveRDS(xiang, file = "data/Xiang_et_al/xiang_fpkm_seurat_object.RDS")
