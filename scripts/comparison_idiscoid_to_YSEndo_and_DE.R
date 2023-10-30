library(Seurat)
library(dplyr)

source("scripts/similarity.R")
source("scripts/utils.R")

###########################################
# 1. Functions for loading markers
###########################################
markerDF2ListHelper <- function(df) {
  markers_list <- as.list(df)
  for (i in 1:length(markers_list)) {
    x <- markers_list[[i]]
    markers_list[[i]] <- convertMouseGeneList(x[x != ""])
  }
  return(markers_list)
}

loadTyserMarkers <- function() {
  #markers <- readRDS(snakemake@input[["Tyser_markers"]])
  markers <- readRDS("results/R_objects/Tyser_markers.RDS")
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markers$gene[markers$cluster == cur_cluster]
    markers_list[[i]] <- cur_markers
  }
  return(markers_list)
}

#######################################################
# 2. Load markers and all genen names for each data set
#######################################################
#Tyser <- readRDS(snakemake@input[["Tyser_seurat_obj"]])
Tyser <- readRDS("results/R_objects/Tyser_seurat_object_processed.RDS")
Tyser_markers <- loadTyserMarkers()
Tyser_all_genes <- row.names(Tyser)

load("data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData")
load("data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData")
load("data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData")
load("data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData")
#load(snakemake@input[["D4_merged"]])
#load(snakemake@input[["D5_20220330_Aug22Revisions"]])
#load(snakemake@input[["mmBmK_D2"]])
#load(snakemake@input[["mmBmK_D3"]])

# Output folder
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision')
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision/figures')
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision/figure_data')
dir.create('results/figures/comparison_idiscoid_to_YSEndo_and_DE')
dir.create('results/figures/comparison_idiscoid_to_YSEndo_and_DE/figures')
dir.create('results/figures/comparison_idiscoid_to_YSEndo_and_DE/figure_data')

####################################################################################
# 3. Perform diff analysis for YSEndo VS DE(P)
####################################################################################
tyser_diff <- FindMarkers(Tyser, ident.1 = "YS Endoderm", ident.2 = "DE(P)")

####################################################################################
# 3. Perform comparisons 
####################################################################################
obj_names <- c("D4_merged","D5_20220330_Aug22Revisions","mmBmK_D2","mmBmK_D3")
display_names <- c("D4_merged","D5","mmBmK_D2","mmBmK_D3")

for(i in 1:length(obj_names)){    
    obj_name <- obj_names[i]
    display_name <- display_names[i]
    
    marker_df <- paste0(obj_name,".markers") %>% as.symbol() %>% eval()
    marker_df <- marker_df[marker_df$p_val_adj < 0.05,]
    res <- list()
    tyser_fc <- list()
    
    # Go through each cluster and compare to DE and YSEndo of Tyser dataset
    jaccard_table <- list()
    for (cluster in unique(marker_df$cluster)){
        cluster_genes <- marker_df$gene[marker_df$cluster == cluster]
        
        # Compute jaccard of unique overlap between cluster VS DE and cluster VS YSEndo.
        DEP_markers <- Tyser_markers[["DE(P)"]]
        YSEndo_markers <- Tyser_markers[["YS Endoderm"]]
        o1 <- setdiff(intersect(DEP_markers, cluster_genes), YSEndo_markers)
        o2 <- setdiff(intersect(YSEndo_markers, cluster_genes), DEP_markers)
        j1 <- length(o1)/length(union(DEP_markers, cluster_genes)) 
        j2 <- length(o2)/length(union(YSEndo_markers, cluster_genes))
        
        # Compute a null distribution of jaccard.
        background_idiscoid <- obj_name %>% as.symbol() %>% eval() %>% rownames()
        background_tyser <- Tyser_all_genes
        size_DEP <- length(DEP_markers)
        size_YSEndo <- length(YSEndo_markers)
        size_idiscoid <- length(cluster_genes)
        distr = c()
        
        # Keep jaccard results for each cluster
        jaccard_table[[length(jaccard_table)+1]] <- c(length(DEP_markers),length(YSEndo_markers),length(cluster_genes),
                                                      length(o1),length(o2),cluster)
        
        for(i in 1:10000){
            DEP_sample <- sample(background_tyser, size_DEP, replace = F)
            YSEndo_sample <- sample(background_tyser, size_YSEndo, replace = F)
            idiscoid_sample <- sample(background_idiscoid, size_idiscoid, replace = F)
            o1_sample <- setdiff(intersect(DEP_sample, idiscoid_sample), YSEndo_sample)
            o2_sample <- setdiff(intersect(YSEndo_sample, idiscoid_sample), DEP_sample)
            j1_sample <- length(o1)/length(union(DEP_sample, idiscoid_sample))
            j2_sample <- length(o2)/length(union(YSEndo_sample, idiscoid_sample))
            distr = c(j2_sample - j1_sample,distr)
        }
        
        # Compute pval of jaccard difference
        pval = length(distr[distr >= j2-j1])/10000
        res[[length(res) + 1]] <- c(cluster, pval)
        
        # Mark the cluster marker genes in tyser_diff dataframe
        tyser_diff$is_marker <- ifelse(rownames(tyser_diff) %in% cluster_genes,"Y","N")
        tyser_diff$cluster <- cluster
        tyser_fc[[length(tyser_fc) + 1]] <- tyser_diff
      }
    
    # Plot jaccard difference
    res_tab <- do.call(rbind,res)
    res_tab <- as.data.frame(res_tab)
    colnames(res_tab) <- c("cluster","pval")
    res_tab$pval <- as.numeric(res_tab$pval)
    res_tab$pval <- p.adjust(res_tab$pval, method = "BH")
    res_tab$pval <- -1 * log10(res_tab$pval + 0.000001)
    
    jaccard_table <- do.call(rbind,jaccard_table)
    colnames(jaccard_table) <- c("Tyser_DEP_markers",
                                "Tyser_YSEndo_markers",
                                "iDisocid_cluster_markers",
                                "overlap(Tyser_DEP_markers, iDisocid_cluster_markers)",
                                "overlap(Tyser_YSEndo_markers, iDisocid_cluster_markers)",
                                "Cluster"
                               )
    
    g1 <- ggplot(res_tab, aes_string("cluster", "pval", fill = "cluster"))
    g1 <- g1 + geom_col()
    g1 <- g1 + geom_hline(yintercept = -1 * log10(0.05))
    g1 <- g1 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))
    g1 <- g1 + labs(x = "iDiscoid clusters", y = "-log10(Pval)", title = "YSEndo - DE(P) similarity ")
    
    # Plot Tyser vocano plot
    tyser_fc <- do.call(rbind,tyser_fc)
    tyser_fc$`-log10(p_val_adj)` <- -1 * log10(tyser_fc$p_val_adj)
    g2 <- ggplot(tyser_fc, aes_string("avg_log2FC", "-log10(p_val_adj)",color = "is_marker"))
    g2 <- g2 + geom_point()
    g2 <- g2 + geom_hline(yintercept = -1 * log10(0.05))
    g2 <- g2 + facet_grid(cols = vars(cluster))
    g2 <- g2 + scale_color_manual(values=c('grey70','red1'))
    g2 <- g2 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none", plot.title = element_text(hjust = 0.5))
    g2 <- g2 + labs(x = "Average log2FC for YSEndo VS DE(P)", y = "-log10(Pval)", title = "YSEndo VS DE(P) genes")
    
    # Save plots, data table and background gene lists.
    outdir = "results/figures/comparison_idiscoid_to_YSEndo_and_DE"
    outname1 <- paste0(outdir,"/figures/",display_name,"_jaccard_comparison.png")
    outname2 <- paste0(outdir,"/figures/",display_name,"_volcano.png")
    outname3 <- paste0(outdir,"/figure_data/",display_name,"_cluster_genes.csv")
    outname4 <- paste0(outdir,"/figure_data/",display_name,"_jaccard_gene_num.csv")
    
    # Save data table
    write.csv(tyser_fc, outname3, row.names = T, col.names = T)
    
    # save jaccard info
    write.table(jaccard_table, outname4, sep = ",")
    
    ggsave(outname1, plot=g1, width = 10*length(unique(marker_df$cluster))/12, height = 5, units = "in") # For better visualization, automatic adjust the width to adapt to number of clusters.

    ggsave(outname2, plot=g2, width = 10*length(unique(marker_df$cluster))/4, height = 5, units = "in") # For better visualization, automatic adjust the width to adapt to number of clusters.
}

file.create(snakemake@output[["outfile"]])