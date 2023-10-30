library(Seurat)
library(magrittr)
library(Matrix)
library(tidyverse)
library(ggplot2)

file_names = unlist(snakemake@input)

# Genes of interest
gene_GATA6 = c("GATA6","GATA4","PDGFRA","SOX17","NODAL","BMP4","HHEX","CER1","LHX1","LEFTY1","LEFTY2")
gene_WT = c("POU5F1","NANOG","SOX2","ISL1","TFAP2A","GATA3","TBXT","MIXL1","GSC","NODAL","BMP4")

# The WT cluster IDs at each time point
WT_cluster_ID <- list("Day 0" = c(0, 1),
                        "Day 0.5" = c(3),
                        "Day 1" = c(3),
                        "Day 1.5" = c(3),
                        "Day 2" = c(4),
                        "Day 3" = c(0),
                        "Day 4" = c(2, 6, 7),
                        "Day 5" = c(0, 3)
                       )

# Load seurat objects
for(fn in file_names){
    load(fn)
}

# Merge all objects
merged <- merge(mmBmK_D0Ri, y = c(mmBmK_12hr, mmBmK_D1, mmBmK_36hr, mmBmK_D2, mmBmK_D3, D4_merged, D5_20220330_Aug22Revisions))

# Change time point names
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D0Ri", "Day 0", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_12hr", "Day 0.5", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D1", "Day 1", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_36hr", "Day 1.5", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D2", "Day 2", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_mmBmK_D3", "Day 3", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_D4_merged", "Day 4", merged@meta.data$orig.ident)
merged@meta.data$orig.ident <- gsub("10X_D5_20220330_Aug22Revisions", "Day 5", merged@meta.data$orig.ident)

# Perform sctransform to remove sequencing depth bias
merged <- SCTransform(merged, verbose = T)

meta <- merged@meta.data
merged@active.assay <- "SCT" # Make sure SCT transformed expressions are used
expr <- GetAssayData(merged, slot = "data")

# Use this to save results
res_GATA6 <- list("Gene"=c(),
            "Average Expression"=c(),
            "Time"=c(),
            "Percent Expressed (FC)"=c()
           )
res_WT <- list("Gene"=c(),
            "Average Expression"=c(),
            "Time"=c(),
            "Percent Expressed (FC)"=c()
           )

# Get mean gene expressions and percentage of gene expressions
for(i in 1:length(WT_cluster_ID)){
    tp <- names(WT_cluster_ID)[i]
    cl <- WT_cluster_ID[[tp]]
    cell_WT <- rownames(meta[(meta$orig.ident == tp) & (meta$seurat_clusters %in% cl),])
    cell_GATA6 <- rownames(meta[(meta$orig.ident == tp) & !(meta$seurat_clusters %in% cl),])
    
    # Some gene might be filtered out in previous steps. These are considered as 'zeros'.
    expr_WT <- expr[intersect(gene_WT,rownames(expr)), cell_WT]
    filtered <- setdiff(gene_WT,rownames(expr))
    expr_zero <- Matrix(0,nrow=length(filtered),ncol = ncol(expr_WT))
    rownames(expr_zero) <- filtered
    expr_WT <- rbind(expr_WT, expr_zero)
    
    expr_GATA6 <- expr[intersect(gene_GATA6,rownames(expr)), cell_GATA6]
    filtered <- setdiff(gene_GATA6,rownames(expr))
    expr_zero <- Matrix(0,nrow=length(filtered),ncol = ncol(expr_GATA6))
    rownames(expr_zero) <- filtered
    expr_GATA6 <- rbind(expr_GATA6, expr_zero)
    
    # compute mean expressions and percentage of gene expressions
    if(ncol(expr_WT) > 0){
        means_WT <- rowMeans(expr_WT)
        percent_WT <- rowSums(expr_WT != 0)/ncol(expr_WT)
    }else{
        means_WT <- rep(0, length(gene_WT))
        percent_WT <- rep(0, length(gene_WT))
        names(means_WT) <- gene_WT
        names(percent_WT) <- gene_WT
    }
    
    if(ncol(expr_GATA6) > 0){
        means_GATA6 <- rowMeans(expr_GATA6)
        percent_GATA6 <- rowSums(expr_GATA6 != 0)/ncol(expr_GATA6)
    }else{
        means_GATA6 <- rep(0, length(gene_GATA6))
        percent_GATA6 <- rep(0, length(gene_GATA6))
        names(means_GATA6) <- gene_GATA6
        names(percent_GATA6) <- gene_GATA6
    }
    
    res_WT$Gene <- c(res_WT$Gene, gene_WT)
    res_WT$`Average Expression` <- c(res_WT$`Average Expression`, means_WT)
    res_WT$Time <- c(res_WT$Time, rep(tp, length(gene_WT)))
    res_WT$`Percent Expressed (FC)` <- c(res_WT$`Percent Expressed (FC)`, percent_WT)
    
    res_GATA6$Gene <- c(res_GATA6$Gene, gene_GATA6)
    res_GATA6$`Average Expression` <- c(res_GATA6$`Average Expression`, means_GATA6)
    res_GATA6$Time <- c(res_GATA6$Time, rep(tp, length(gene_GATA6)))
    res_GATA6$`Percent Expressed (FC)` <- c(res_GATA6$`Percent Expressed (FC)`, percent_GATA6)
}

res_WT = as_tibble(res_WT)
res_GATA6 = as_tibble(res_GATA6)

# Transform percent expressed by each gene's mean percent expressed in all time points
for(gene in unique(res_WT$Gene)){
    mean_percent =  mean(res_WT[res_WT$Gene == gene,]$`Percent Expressed (FC)`)
    if(mean_percent != 0){
        res_WT[res_WT$Gene == gene,]$`Percent Expressed (FC)` = res_WT[res_WT$Gene == gene,]$`Percent Expressed (FC)`/mean_percent
    }
}

for(gene in unique(res_GATA6$Gene)){
    mean_percent =  mean(res_GATA6[res_GATA6$Gene == gene,]$`Percent Expressed (FC)`)
    if(mean_percent != 0){
        res_GATA6[res_GATA6$Gene == gene,]$`Percent Expressed (FC)` = res_GATA6[res_GATA6$Gene == gene,]$`Percent Expressed (FC)`/mean_percent
    }
}

# Create output folder
dir.create('results/figures/expression_change')

# Generate expression change plots
plot_theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 15, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 15)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 15)),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  theme(plot.title = element_text(hjust = 0.5)),
  scale_size_continuous(range = c(0, 15)),
  scale_x_discrete(position = "bottom")
)

WT <- ggplot(data = res_WT,
       aes(Gene, Time)) + geom_point(aes(size = `Percent Expressed (FC)`,color = `Average Expression`)
       ) + ggtitle("WT Clusters") + plot_theme

GATA6 <- ggplot(data = res_GATA6,
       aes(Gene, Time)) + geom_point(aes(size = `Percent Expressed (FC)`,color = `Average Expression`)
       ) + ggtitle("iGATA6 Clusters") + plot_theme

ggsave(snakemake@output[["WT"]], WT, width = 12, height = 7)
ggsave(snakemake@output[["GATA6"]],GATA6,width = 12, height = 7)