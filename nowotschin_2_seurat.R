library(Seurat)
library(data.table)

# counts <- read.csv("nowotschin_data/sc_endoderm_all_cells_counts.csv", header = TRUE, row.names = 1)
counts <- fread(file = "nowotschin_data/sc_endoderm_all_cells_counts.csv")
meta <- read.csv("nowotschin_data/sc_endoderm_all_cells_metadata.csv", header = TRUE, row.names = 1)

saveRDS(counts, file = "nowotschin_data/sc_endoderm_all_cells_counts.RDS")
saveRDS(meta, file = "nowotschin_data/sc_endoderm_all_cells_metadata.RDS")

# Get E6.5
counts <- counts[meta$Timepoint == "E6.5",]
setDF(counts)
tmp <- counts[, -1]
rownames(tmp) <- counts[, 1]
counts <- tmp
meta <- meta[meta$Timepoint == "E6.5",]
nowotschin_e6.5 <- CreateSeuratObject(counts = t(counts), project = "NowotschinE6.5", assay = "RNA", meta.data = meta)
nowotschin_e6.5.filtered <- CreateSeuratObject(counts = t(counts), project = "NowotschinE6.5", assay = "RNA", meta.data = meta, min.cells = 5, min.features = 200)

saveRDS(nowotschin_e6.5, file = "nowotschin_data/nowotschin_E6.5_seurat_object.RDS")
saveRDS(nowotschin_e6.5.filtered, file = "nowotschin_data/nowotschin_E6.5_filtered_seurat_object.RDS")
