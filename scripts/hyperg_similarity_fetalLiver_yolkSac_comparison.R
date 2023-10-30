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

markersDF2List <- function(markersDF) {
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


#######################################################
# 2. Load markers and all genen names for each data set
#######################################################

feloDeslo_seurat_obj_names <- c("DesLO_D17",
                             "coculture_FeLO",
                             "FeLO_monoculture_cytokine_minus",
                             "FeLO_monoculture_cytokine_plus",
                             "FeLO_D17_JH"
                            )

compared_seurat_obj_names <- c("AFLYS",
                              "yolk_sac",
                              "YSFetalLiverFullMetadata",
                              "ZhaiCynomolgus"
                              )

# Load all data (RData)
for(i in 1:7){
    load(snakemake@input[[i]])
}

# Load all data (rds)
DesLO_D17 <- readRDS(snakemake@input[[8]])
FeLO_D17_JH <- readRDS(snakemake@input[[9]])
YSFetalLiverFullMetadata.markers <- readRDS(snakemake@input[[10]])
DesLO_D17.markers <- readRDS(snakemake@input[[11]])
FeLO_D17_JH.markers <- readRDS(snakemake@input[[12]])

# Output folder
#dir.create('results/figures/hypergeometric_comparisons/fetalLiver_yolksac')
#dir.create('results/figures/hypergeometric_comparisons/fetalLiver_yolksac/figures')
#dir.create('results/figures/hypergeometric_comparisons/fetalLiver_yolksac/figure_data')
dir.create('results/figures/hypergeometric_comparisons/figure_data')
table <- list()
####################################################################################
# 3. Perform comparisons. Generate background genes and perform hypergeo comparison 
####################################################################################
for(i in 1:length(feloDeslo_seurat_obj_names)){
    feloDeslo_obj_name <- feloDeslo_seurat_obj_names[i]
    
    feloDeslo_seurat_obj <- eval(as.symbol(feloDeslo_obj_name))
    feloDeslo_gene_names <- rownames(feloDeslo_seurat_obj)
    feloDeslo_markers_df <- paste0(feloDeslo_obj_name,".markers") %>% as.symbol() %>% eval()
    feloDeslo_markers_df <- feloDeslo_markers_df[feloDeslo_markers_df$p_val_adj < 0.05,]
    feloDeslo_cluster_names <- unique(feloDeslo_markers_df$cluster)
    
    # Get markers in each cluster for feloDeslo objects
    feloDeslo_markers <- markersDF2List(feloDeslo_markers_df)
    
    # Compare to other data sets
    for(k in 1:length(compared_seurat_obj_names)){
        compared_obj_name = compared_seurat_obj_names[k]
        
        #ZhaiCynomolgusEmbryoProcessed, YolkSacFetalLiver, and AdultPlusFetalLiverYolkSac were split by days and then comapred to feloDeslo
        compared_seurat_obj <- eval(as.symbol(compared_obj_name))
        compared_gene_names <- row.names(compared_seurat_obj)
        
        if(compared_obj_name == "yolk_sac"){
            compared_markers_df <- paste0(compared_obj_name,".markers") %>% as.symbol %>% eval()
            compared_markers <- markersDF2List(compared_markers_df)
            background <- intersect(feloDeslo_gene_names, compared_gene_names)
            
            # The actual comparison by hypergeometric test
            scores <- compareMarkersBetweenDatasets(
                                   feloDeslo_markers,
                                   compared_markers,
                                   background,
                                   feloDeslo_obj_name,
                                   compared_obj_name
            )
            #g + labs(x = paste(compared_obj_name,"Cell Type"), y = "-log10(Hyperg Pval)")
            
            # Save plot and plot data
            #outdir = "results/figures/hypergeometric_comparisons/fetalLiver_yolksac/"
            #outname <- paste0(outdir,
            #                   "figures/",
            #                   feloDeslo_obj_name,"-VS-",
            #                   compared_obj_name,
            #                   ".png"
            #                  )

            #ggsave(outname, width = 10*length(feloDeslo_markers)/4, height = 5, units = "in") # For better visualization, automatic adjust the width so that labels won't overlap with each other.
            #outname <- paste0(outdir,
            #                   "figure_data/",
            #                   feloDeslo_obj_name,"-VS-",
            #                   compared_obj_name,
            #                   ".csv"
            #
            outdir <- "results/figures/hypergeometric_comparisons/figure_data/"
            outname <- paste0(outdir,
                               feloDeslo_obj_name,"-VS-",
                               compared_obj_name,
                               ".csv"
                              )
                                
            write.csv(scores, outname, row.names = F)           
            #write.csv(ggplot_build(g)$plot$data, outname, row.names = F)
            
        # For these datasets, we compared our data to them in different days separately
        }else{
            
            if(compared_obj_name == "AFLYS"){
                dataset_name = "AFLYS"
                day_col = "age"
            }else if(compared_obj_name == "YSFetalLiverFullMetadata"){
                dataset_name = "YSEL"
                day_col = "stage"
            }else if(compared_obj_name == "ZhaiCynomolgus"){
                dataset_name = "ZhaiCynomolgus"
                day_col = "stage"
            }
                for(day in unique(compared_seurat_obj@meta.data[[day_col]])){
                    
                    # Get markers for each day
                    day_name <- paste(dataset_name,day,sep = "_")
                    if(exists(day_name %>% paste0(".markers"))){
                        day_markers_df <- day_name %>% paste0(".markers") %>% as.symbol() %>% eval()
                        day_markers_df <- day_markers_df[day_markers_df$p_val_adj < 0.05,]
                        day_markers <- markersDF2List(day_markers_df)
                        background <- intersect(feloDeslo_gene_names, compared_gene_names)

                        # The actual comparison by hypergeometric test
                        #g <- compareMarkersBetweenDatasets(
                        #                       feloDeslo_markers,
                        #                       day_markers,
                        #                       background,
                        #                       feloDeslo_obj_name,
                        #                       day_name
                        #)
                        
                        #g + labs(x = paste(day_name,"Cell Type"), y = "-log10(Hyperg Pval)")
                            
                        # Save plot 
                        #outdir = "results/figures/hypergeometric_comparisons/fetalLiver_yolksac/"
                        #outname <- paste0(outdir,
                        #                   "figures/",
                        #                   feloDeslo_obj_name,"-VS-",
                        #                   day_name,
                        #                   ".png"
                        #                  )

                        #ggsave(outname, width = 10*length(feloDeslo_markers)/4, height = 5, units = "in") # For better visualization, automatic adjust the width so that labels won't overlap with each other.
                        scores <- compareMarkersBetweenDatasets(
                                               feloDeslo_markers,
                                               day_markers,
                                               background,
                                               feloDeslo_obj_name,
                                               day_name
                        )
                        outdir <- "results/figures/hypergeometric_comparisons/figure_data/"
                        outname <- paste0(outdir,
                                           feloDeslo_obj_name,"-VS-",
                                           day_name,
                                           ".csv"
                                          )
                        
                       
                        write.csv(scores, outname, row.names = F)
                    }
                 }
                    
        }
    }
}

# If successful
file.create(snakemake@output[[1]])