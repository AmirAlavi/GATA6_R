# Libraries
library(dplyr)
library(Seurat)

day_markers_list <- list()
gene_sets <- list()
days <- list()
j <- 1
for (data in snakemake@input) {
    day <- basename(data)
    day <- sub("mmBmK", "", day)
    day <- gsub("^_", "", day)
    day <- unlist(strsplit(day, "_"))[[1]]
    print(day)

    load(data)
    day_markers <- get(paste0("mmBmK_", day, ".markers"))
    day_markers <- subset(day_markers, p_val_adj < 0.05)
    print(str(day_markers))
    if (nrow(day_markers) > 0) {
        day_markers_list[[j]] <- day_markers
        days[[j]] <- day
        gene_sets[[j]] <- row.names(get(paste0("mmBmK_", day)))
        j <- j + 1
    }
}
names(day_markers_list) <- days
names(gene_sets) <- days
print(str(day_markers_list))
print(str(gene_sets))
saveRDS(day_markers_list, file = snakemake@output[[1]])
saveRDS(gene_sets, file = snakemake@output[[2]])
