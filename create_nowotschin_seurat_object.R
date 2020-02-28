library(Seurat)
library(data.table)


# This is a large file, so using data.table instead of data.frame
counts <- fread(file = "data/Nowotschin_et_al/sc_endoderm_all_cells_counts.csv")
meta <- read.csv("data/Nowotschin_et_al/sc_endoderm_all_cells_metadata.csv", header = TRUE, row.names = 1)

saveRDS(counts, file = "data/Nowotschin_et_al/sc_endoderm_all_cells_counts.RDS")
saveRDS(meta, file = "data/Nowotschin_et_al/sc_endoderm_all_cells_metadata.RDS")

# Get E6.5 cells
counts <- counts[meta$Timepoint == "E6.5",]
meta <- meta[meta$Timepoint == "E6.5",]
# convert to data.frame
setDF(counts)

# These lines are to mend the rownames (gene names)
tmp_rownames <- counts[, 1]
counts <- counts[, -1]
rownames(counts) <- tmp_rownames

counts <- t(counts)

saveRDS(counts, file = "data/Nowotschin_et_al/nowotschin_E6.5_counts.RDS")
saveRDS(meta, file = "data/Nowotschin_et_al/nowotschin_E6.5_metadata.RDS")

nowotschin_e6.5 <- CreateSeuratObject(counts = counts, project = "NowotschinE6.5", assay = "RNA", meta.data = meta, min.cells = 5, min.features = 200)

saveRDS(nowotschin_e6.5, file = "data/Nowotschin_et_al/nowotschin_E6.5_seurat_object.RDS")
