library(Seurat)
library(data.table)


# This is a large file, so using data.table instead of data.frame
counts <- fread(file = snakemake@input[[1]])
counts <- as.matrix(counts, rownames=1)
meta <- read.csv(snakemake@input[[2]], header = TRUE, row.names = 1)

xiang <- CreateSeuratObject(counts = counts, project = "Xiang", assay = "RNA", meta.data = meta)

saveRDS(xiang, file = snakemake@output[[1]])
