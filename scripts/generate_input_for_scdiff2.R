library(Seurat)
library(Matrix)
library(magrittr)
file_names = unlist(snakemake@input)

# Load seurat objects
for(fn in file_names[1:8]){
    load(fn)
}

# Load cluster labels
cluster_label <- read.table(file_names[9],row.names = 1,sep = ",", header = T)

# The WT cluster IDs at each time point
WT_cluster_ID <- list("0" = c(0, 1),
                        "0.5" = c(3),
                        "1" = c(3),
                        "1.5" = c(3),
                        "2" = c(4),
                        "3" = c(0),
                        "4" = c(2, 6, 7),
                        "5" = c(0, 3)
                       )


# Merge all objects
merged <- merge(mmBmK_D0Ri, y = c(mmBmK_12hr, mmBmK_D1, mmBmK_36hr, mmBmK_D2, mmBmK_D3, D4_merged, D5_20220330_Aug22Revisions))

# Change time point names
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D0Ri", "0", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_12hr", "0.5", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D1", "1", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_36hr", "1.5", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D2", "2", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D3", "3", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_D4_merged", "4", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_D5_20220330_Aug22Revisions", "5", merged@meta.data$orig.ident)

meta <- merged@meta.data
meta$cluster_label <- NA

# At each time point, get cluster labels
for(tp in unique(meta$orig.ident)){
    
    # Day 3 has subclusters that further divides cluster 0
    if(tp == "3"){
        tp_cls <- merged@active.ident[meta$orig.ident == tp] %>% unique %>% as.character
    }else{
        tp_cls <- meta$seurat_clusters[meta$orig.ident == tp] %>% unique
    }
    for(i in 1:ncol(cluster_label)){
        if(!is.na(cluster_label[tp,i]) && cluster_label[tp,i] != ""){
            if(length(grep(",",cluster_label[tp,i]))){
                lab_cls <- cluster_label[tp,i]  %>% strsplit(",") %>% unlist
            }else{
                lab_cls <- cluster_label[tp,i]
            }
            lab <- colnames(cluster_label)[i]
            for(cl in intersect(tp_cls,lab_cls)){
                if(tp == "3"){
                    meta[(meta$orig.ident == tp) & (merged@active.ident == cl),]$cluster_label = lab
                }else{
                    meta[(meta$orig.ident == tp) & (meta$seurat_clusters == cl),]$cluster_label = lab
                }
                
            }
        } 
    }  
}

merged@meta.data <- meta

# At each time point, separate WT clusters and iGATA6 clusters
cell_WT <- c()
cell_GATA6 <- c()
for(i in 1:length(WT_cluster_ID)){
    tp <- names(WT_cluster_ID)[i]
    cl <- WT_cluster_ID[[tp]]
    cell_WT <- c(rownames(meta[(meta$orig.ident == tp) & (meta$seurat_clusters %in% cl),]),cell_WT)
    cell_GATA6 <- c(rownames(meta[(meta$orig.ident == tp) & !(meta$seurat_clusters %in% cl),]),cell_GATA6)
}

merged_WT <- merged[,cell_WT]
merged_GATA6 <- merged[,cell_GATA6]

# Perform sctransform to remove sequencing depth bias
#merged_WT <- SCTransform(merged_WT, verbose = T)
#merged_GATA6 <- SCTransform(merged_GATA6, verbose = T)

# Identify HVG and keep only the HVGs
#FindVariableFeatures(merged_WT, nfeatures = 3000)
#FindVariableFeatures(merged_GATA6, nfeatures = 3000)
#merged_WT <- merged_WT[VariableFeatures(merged_WT),]
#merged_GATA6 <- merged_GATA6[VariableFeatures(merged_GATA6),]

meta_WT <- merged_WT@meta.data
#merged_WT@active.assay <- "SCT" # Make sure SCT transformed expressions are used
meta_GATA6 <- merged_GATA6@meta.data
#merged_GATA6@active.assay <- "SCT" # Make sure SCT transformed expressions are used

# Extract sctransformed expression data
#expr_WT <- GetAssayData(merged_WT, slot = "data")
#expr_GATA6 <- GetAssayData(merged_GATA6, slot = "data")
expr_WT <- GetAssayData(merged_WT, slot = "counts")
expr_GATA6 <- GetAssayData(merged_GATA6, slot = "counts")

# Extract meta data information (need cell ids, time points, cell labels, cluster labels)
cell_id_WT <- rownames(merged_WT@meta.data)
time_point_WT <- merged_WT@meta.data$orig.ident
cell_label_WT <- merged_WT@meta.data$cluster_label
cluster_WT <- merged_WT@meta.data$seurat_clusters

# For day 3 wild type, we need subcluster labels
cluster_WT[merged_WT@meta.data$orig.ident == "3"] <- merged_WT@active.ident[merged_WT@meta.data$orig.ident == "3"] %>% as.character

cell_id_GATA6 <- rownames(merged_GATA6@meta.data)
time_point_GATA6 <- merged_GATA6@meta.data$orig.ident
cell_label_GATA6 <- merged_GATA6@meta.data$cluster_label
cluster_GATA6 <- merged_GATA6@meta.data$seurat_clusters

# Generate input file for scdiff2
gene_names_WT <- rownames(merged_WT)
expr_tab_WT <- cbind(cell_id_WT, time_point_WT,cell_label_WT,cluster_WT,as.data.frame(t(expr_WT)))
gene_names_GATA6 <- rownames(merged_GATA6)
expr_tab_GATA6 <- cbind(cell_id_GATA6, time_point_GATA6,cell_label_GATA6,cluster_GATA6,as.data.frame(t(expr_GATA6)))

dir.create("results/scdiff2")
write.table(expr_tab_WT, snakemake@output[["scdiff2_out_WT"]], row.names = F, sep = "\t", quote = F)
write.table(expr_tab_GATA6, snakemake@output[["scdiff2_out_GATA6"]], row.names = F, sep = "\t", quote = F)