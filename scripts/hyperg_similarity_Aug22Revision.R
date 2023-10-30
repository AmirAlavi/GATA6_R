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

markersDF2List_Xiang <- function(markersDF) {
  n_clusters <- length(unique(markersDF$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markersDF$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markersDF$gene[markersDF$cluster == cur_cluster]
    markers_list[[i]] <- cur_markers
  }
  return(markers_list)
}

loadTyserMarkers <- function() {
  markers <- readRDS(snakemake@input[["Tyser_markers"]])
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

loadMaMarkers <- function() {
  markers <- read.csv(snakemake@input[["Ma_markers"]], colClasses = "character", row.names = 1)
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    markers_list[[i]] <- markers$gene[markers$cluster == cur_cluster]
  }
  return(markers_list)
}

loadNowotschinE4.5Markers <- function() {
  markers <- read.csv(snakemake@input[["Nowotschin_E4pt5_markers"]], colClasses = "character")
  markerDF2ListHelper(markers)
}

loadNowotschinE5.5Markers <- function() {
  markers <- read.csv(snakemake@input[["Nowotschin_E5pt5_markers"]], colClasses = "character")
  markerDF2ListHelper(markers)
}

loadNowotschinE6.5Markers <- function() {
  markers <- readRDS(snakemake@input[["Nowotschin_E6pt5_markers"]])
  n_clusters <- length(unique(markers$cluster))
  markers_list <- vector(mode = "list", length = n_clusters)
  names(markers_list) <- unique(markers$cluster)
  for (i in 1:length(markers_list)) {
    cur_cluster <- names(markers_list)[i]
    cur_markers <- markers$gene[markers$cluster == cur_cluster]
    markers_list[[i]] <- convertMouseGeneList(cur_markers)
  }
  return(markers_list)
}

loadNowotschinE7.5Markers <- function() {
  markers <- read.csv(snakemake@input[["Nowotschin_E7pt5_markers"]], colClasses = "character")
  markerDF2ListHelper(markers)
}

loadSalaMarkers <- function() {
  markers <- read.csv(snakemake@input[["Sala_markers"]], colClasses = "character")
  markerDF2ListHelper(markers)
}

#######################################################
# 2. Load markers and all genen names for each data set
#######################################################
Tyser <- readRDS(snakemake@input[["Tyser_seurat_obj"]])
Tyser_markers <- loadTyserMarkers()
Tyser_all_genes <- row.names(Tyser)
rm(Tyser)

Ma_markers <- loadMaMarkers()
Ma_counts <- read.csv(snakemake@input[["Ma_counts"]], row.names = 1)
Ma_all_genes <- rownames(Ma_counts)

Nowotschin_E4.5_markers <- loadNowotschinE4.5Markers()
Nowotschin_E5.5_markers <- loadNowotschinE5.5Markers()
Nowotschin_E6.5_markers <- loadNowotschinE6.5Markers()
Nowotschin_E7.5_markers <- loadNowotschinE7.5Markers()
Nowotschin <- readRDS(snakemake@input[["Nowotschin_E6pt5_seurat_obj"]])
Nowotschin_all_genes <- row.names(Nowotschin)
Nowotschin_all_genes <- convertMouseGeneList(Nowotschin_all_genes)
Nowotschin_E4.5_all_genes <- Nowotschin_all_genes
Nowotschin_E5.5_all_genes <- Nowotschin_all_genes
Nowotschin_E6.5_all_genes <- Nowotschin_all_genes
Nowotschin_E7.5_all_genes <- Nowotschin_all_genes
rm(Nowotschin)

Xiang <- readRDS(snakemake@input[["Xiang_seurat_obj"]])
Xiang_markers <- readRDS(snakemake@input[["Xiang_markers"]])
Xiang_all_genes <- row.names(Xiang)
rm(Xiang)

Sala_markers <- loadSalaMarkers()
Sala_all_genes <- read.csv(snakemake@input[["Sala_all_genes"]], header = FALSE, sep = "\t", colClasses = "character", col.names = c("Ensembl", "Symbol"))$Symbol ##############input
Sala_all_genes <- convertMouseGeneList(Sala_all_genes)

# Load d0-d7 data
for(i in 14:25){
    load(snakemake@input[[i]])
}

# Load aug22 data
for(i in 26:33){
    load(snakemake@input[[i]])
}

other_dataset_names <- c("Tyser",
                          "Ma",
                          "Nowotschin_E4pt5",
                          "Nowotschin_E5pt5",
                          "Nowotschin_E6pt5",
                          "Nowotschin_E7pt5",
                          "Xiang",
                          "Sala"
                         ) 
other_dataset_display_names <- c(
                        "Tyser",
                          "Ma",
                          "Nowotschin_E4.5",
                          "Nowotschin_E5.5",
                          "Nowotschin_E6.5",
                          "Nowotschin_E7.5",
                          "Xiang",
                          "Sala"
                        )

aug22_seurat_obj_names <- names(snakemake@input)[59:66]
aug22_display_names <- c("D4_merged",
                         "D4_SB431542",
                         "D5",
                         "D7",
                         "FeLO_monoculture_cytokine_minus",
                         "FeLO_monoculture_cytokine_plus",
                         "D8_H1iDiscoid",
                         "yolk_sac"
                        )

D0ToD8_seurat_obj_names <- c("mmBmK_12hr",
                             "mmBmK_36hr",
                             "mmBmK_3hr",
                             "mmBmK_D0D1",
                             "mmBmK_D0D5",
                             "mmBmK_D0Ri",
                             "mmBmK_D1",
                             "mmBmK_D2",
                             "mmBmK_D3",
                             "D4_merged",
                             "D4_SB431542_Aug22Revisions",
                             "D5_20220330_Aug22Revisions",
                             "D7_20220307_Aug22Revisions",
                             "D8_H1iDiscoid_Aug22Revision"
                            )

D0ToD8_display_names <- c("mmBmK_12hr",
                          "mmBmK_36hr",
                          "mmBmK_3hr",
                          "mmBmK_D0D1",
                          "mmBmK_D0D5",
                          "mmBmK_D0Ri",
                          "mmBmK_D1",
                          "mmBmK_D2",
                          "mmBmK_D3",
                          "D4_merged",
                          "D4_SB431542",
                          "D5",
                          "D7",
                          "D8_H1iDiscoid"
                          )

D22_seurat_obj_names <- c("FeLO_monoculture_cytokine_minus","FeLO_monoculture_cytokine_plus")

# Output folder
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision')
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision/figures')
#dir.create('results/figures/hypergeometric_comparisons/aug22Revision/figure_data')
dir.create('results/figures/hypergeometric_comparisons/figure_data')
####################################################################################
# 3. Perform comparisons. Generate background genes and perform hypergeo comparison 
####################################################################################
for(i in 1:length(aug22_seurat_obj_names)){
    
    aug22_obj_name <- aug22_seurat_obj_names[i]
    aug22_display_name <- aug22_display_names[i]
    
    # Get markers and clusters for aug22 revised data
    aug22_seurat_obj <- eval(as.symbol(aug22_obj_name))
    aug22_gene_names <- row.names(aug22_seurat_obj)
    aug22_markers_df <- paste0(aug22_obj_name,".markers") %>% as.symbol() %>% eval()
    aug22_markers_df <- aug22_markers_df[aug22_markers_df$p_val_adj < 0.05,]
    aug22_cluster_names <- unique(aug22_markers_df$cluster)
    
    # Get markers in each cluster for aug22 revised data
    aug22_markers <- vector(mode = "list", length = length(aug22_cluster_names))
    names(aug22_markers) <- aug22_cluster_names
    for (j in 1:length(aug22_markers)) {
      cluster <- names(aug22_markers)[j]
      aug22_markers[[j]] <- aug22_markers_df$gene[aug22_markers_df$cluster == cluster]
    }
    
    # Compare to embryo data sets
    for(k in 1:length(other_dataset_names)){
        
        other_dataset_name = other_dataset_names[k]
        other_dataset_display_name = other_dataset_display_names[k]
        
        # Markers in Xiang data set are split into different days
        if(other_dataset_name == "Xiang"){
            
            Xiang_gene_names <- paste0(other_dataset_name,"_all_genes") %>% as.symbol() %>% eval()
            Xiang_markers <- paste0(other_dataset_name,"_markers") %>% as.symbol() %>% eval()
            for(x in 1:length(Xiang_markers)){
                day <- names(Xiang_markers)[x]
                day_markers <- markersDF2List_Xiang(Xiang_markers[[day]])
                background <- intersect(aug22_gene_names, Xiang_gene_names)
                scores <- compareMarkersBetweenDatasets(aug22_markers,
                                           day_markers,
                                           background,
                                           aug22_display_name,
                                           paste0("Xiang_",gsub(" ","_",day,perl = T),"_cluster"))
                
                #g + labs(x = paste("Xiang",day,"Cell Type"), y = "-log10(Hyperg Pval)")
                
                # Save plot
                #outdir = "results/figures/hypergeometric_comparisons/figure_data/"
                #outname <- paste0(outdir,
                #                   aug22_display_name,"-VS-",
                #                   "Xiang_",
                #                   day,
                #                   ".png"
                #                  )
                #ggsave(outname, width = 10*length(aug22_markers)/4, height = 5, units = "in") # For better visualization, automatic adjust the width so that labels won't overlap with each other.
                # save tables
                outdir = "results/figures/hypergeometric_comparisons/figure_data/"
                outname <- paste0(outdir,
                                   aug22_display_name,"-VS-",
                                   "Xiang_",
                                  day,
                                   ".csv"
                                  )
                write.csv(scores, outname, row.names = F)
            }
        }else{
            
            other_dataset_name <- ifelse(any(grep("Nowotschin",other_dataset_name)),sub("pt",".",other_dataset_name),other_dataset_name)
            other_gene_names <- paste0(other_dataset_name,"_all_genes") %>% as.symbol() %>% eval()  
            other_markers <- paste0(other_dataset_name,"_markers") %>% as.symbol() %>% eval()
            background <- intersect(aug22_gene_names, other_gene_names)
            
            scores <- compareMarkersBetweenDatasets(aug22_markers,
                                               other_markers,
                                               background,
                                               aug22_display_name,
                                               paste0(other_dataset_display_name, "_cluster"))
            #g + labs(x = paste0(other_dataset_display_name, " Cell Type"), y = "-log10(Hyperg Pval)")
            
            # Save plot
            #outdir = "results/figures/hypergeometric_comparisons/aug22Revision/"
            #outname <- paste0(outdir,
            #                   "figures/",
            #                   aug22_display_name,"_VS_",
            #                   other_dataset_display_name,
            #                   ".png"
            #                  )
            #ggsave(outname, width = 10*length(aug22_markers)/4, height = 5, units = "in")
            
            # save table
            outname <- paste0("results/figures/hypergeometric_comparisons/figure_data/",
                               aug22_display_name,"-VS-",
                               other_dataset_display_name,
                               ".csv"
                              )
            write.csv(scores, outname, row.names = F)
        }
    }
    
    # Compare to D0-D8 datasets
    for(k in 1:length(D0ToD8_seurat_obj_names)){
        
        D0ToD8_obj_name <- D0ToD8_seurat_obj_names[k]
        D0ToD8_display_name <- D0ToD8_display_names[k]
        
        # Get markers and clusters
        D0ToD8_seurat_obj <- eval(as.symbol(D0ToD8_obj_name))
        D0ToD8_gene_names <- as.symbol(D0ToD8_obj_name) %>% eval() %>% row.names()
        D0ToD8_markers_df <- paste0(D0ToD8_obj_name,".markers") %>% as.symbol() %>% eval()
        D0ToD8_markers_df <- D0ToD8_markers_df[D0ToD8_markers_df$p_val_adj < 0.05,]
        D0ToD8_cluster_names <- unique(D0ToD8_markers_df$cluster) 
        
        # Get markers in each cluster
        D0ToD8_markers <- vector(mode = "list", length = length(D0ToD8_cluster_names))
        names(D0ToD8_markers) <- D0ToD8_cluster_names
        for (x in 1:length(D0ToD8_markers)) {
          cluster <- names(D0ToD8_markers)[x]
          D0ToD8_markers[[x]] <- D0ToD8_markers_df$gene[D0ToD8_markers_df$cluster == cluster]
        }
        
        background <- intersect(D0ToD8_gene_names, aug22_gene_names)
        scores <- compareMarkersBetweenDatasets(aug22_markers,
                                           D0ToD8_markers,
                                           background,
                                           aug22_display_name,
                                           paste0(D0ToD8_display_name, "_cluster"))
        #g + labs(x = paste0(D0ToD8_display_name, " Cell Type"), y = "-log10(Hyperg Pval)")
        
        # Save plot
        #outdir = "results/figures/hypergeometric_comparisons/aug22Revision/"
        #outname <- paste0(outdir,
        #                   "figures/",
        #                   aug22_display_name,"_VS_",
        #                   D0ToD8_display_name,
        #                   ".png"
        #                  )
        #ggsave(outname, width = 10*length(aug22_markers)/4, height = 5, units = "in")
        
        # save table
        outdir = "results/figures/hypergeometric_comparisons/figure_data/"
        outname <- paste0(outdir,
                           aug22_display_name,"-VS-",
                           D0ToD8_display_name,
                           ".csv"
                          )
        write.csv(scores, outname, row.names = F)
    }
    
    # Additonally, compare yolk sac data set to D22 data sets
    if(aug22_obj_name == "yolk_sac"){
        for(D22_seurat_obj_name in D22_seurat_obj_names){
            D22_seurat_obj <- eval(as.symbol(D22_seurat_obj_name))
            D22_gene_names <- row.names(D22_seurat_obj)
            D22_markers_df <- paste0(D22_seurat_obj_name,".markers") %>% as.symbol() %>% eval()
            D22_markers_df <- D22_markers_df[D22_markers_df$p_val_adj < 0.05,]
            D22_cluster_names <- unique(D22_markers_df$cluster)

            # Get markers in each cluster for aug22 revised data
            D22_markers <- vector(mode = "list", length = length(D22_cluster_names))
            names(D22_markers) <- D22_cluster_names
            
            for (k in 1:length(D22_markers)) {
              cluster <- names(D22_markers)[k]
              D22_markers[[k]] <- D22_markers_df$gene[D22_markers_df$cluster == cluster]
            }
            background <- intersect(aug22_gene_names, D22_gene_names)
            
            scores <- compareMarkersBetweenDatasets(aug22_markers,
                                               D22_markers,
                                               background,
                                               aug22_display_name,
                                               paste0(D22_seurat_obj_name, "_cluster"))
            #g + labs(x = paste0(D22_seurat_obj_name, " Cell Type"), y = "-log10(Hyperg Pval)")
            
            # Save plot
            #outdir = "results/figures/hypergeometric_comparisons/aug22Revision/"
            #outname <- paste0(outdir,
            #                   "figures/",
            #                   aug22_display_name,"_VS_",
            #                   D22_seurat_obj_name,
            #                   ".png"
            #                  )
            #ggsave(outname, width = 10*length(aug22_markers)/4, height = 5, units = "in")
            
            # save table
            outdir = "results/figures/hypergeometric_comparisons/figure_data/"
            outname <- paste0(outdir,
                               aug22_display_name,"_VS_",
                               D22_seurat_obj_name,
                               ".csv"
                              )
            write.csv(scores, outname, row.names = F)
        }
    }
}
file.create(snakemake@output[["outfile"]])