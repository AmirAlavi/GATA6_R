# Libraries
library(dplyr)
library(Seurat)

xiang <- readRDS(snakemake@input[[1]])
day_markers_list <- list()
all_days <- unique(xiang$Age)
days <- list()
j <- 1
for (day in all_days) {
    print(day)
    day_data <- subset(x = xiang, subset = Age == day)
    print(dim(day_data))
    day_markers <- FindAllMarkers(object = day_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    day_markers <- subset(day_markers, p_val_adj < 0.05)
    if (nrow(day_markers) > 0) {
        day_markers_list[[j]] <- day_markers
        days[[j]] <- day
        j <- j + 1
    }
}
names(day_markers_list) <- days
print(str(day_markers_list))
saveRDS(day_markers_list, file = snakemake@output[[1]])
