library(ggplot2)
source("similarity.R")

# https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x, verbose = FALSE){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])

  if (verbose) {
    # Print the first 6 genes found to the screen
    print(head(humanx))
  }
  return(humanx)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x, verbose = FALSE){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  if (verbose) {
    # Print the first 6 genes found to the screen
    print(head(humanx))
  }
  return(humanx)
}

# analyze the similarity of the within-dataset cluster markers
compareMarkersWithinDataset <- function(markers) {
  n_scores <- 2 * length(markers)^2
  jaccard_sims <- numeric(n_scores)
  clusters1 <- numeric(n_scores)
  clusters2 <- numeric(n_scores)
  part <- character(n_scores)

  n <- 1
  for(i in 1:length(markers)) {
    cluster1 <- names(markers)[i]
    print(cluster1)
    cluster1_markers <- markers[[i]]
    for(j in 1:length(markers)) {
      cluster2 <- names(markers)[j]
      print(cluster2)
      cluster2_markers <- markers[[j]]
      
      jaccard <- similarityJaccard(cluster1_markers, cluster2_markers, verbose = TRUE)
      
      clusters1[n] <- cluster1
      clusters2[n] <- cluster2
      jaccard_sims[n] <- jaccard
      part[n] <- "similarity"
      
      clusters1[n+1] <- cluster1
      clusters2[n+1] <- cluster2
      jaccard_sims[n+1] <- 1 - jaccard
      part[n+1] <- "remainder"
      
      n <- n + 2
    }
  }
  
  scores <- data.frame("A" = clusters1, "B" = clusters2, "Jaccard_Similarity" = jaccard_sims, "Partition" = part)
  g <- ggplot(scores, aes(x=1, y=Jaccard_Similarity, fill = Partition))
  g <- g + geom_col(color = 'black')
  g <- g + coord_polar("y", start=0)
  g <- g + facet_grid(rows = vars(B), cols = vars(A))
  g + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
}

compareMarkersBetweenDatasets <- function(markers1, markers2, background, name_markers1, name_markers2, species1 = c("human", "mouse"), species2 = c("human", "mouse")) {
  n_scores = length(markers1) * length(markers2)
  scores <- numeric(n_scores)
  clusters1 <- character(n_scores)
  clusters2 <- character(n_scores)
  
  n <- 1
  for(i in 1:length(markers1)){
    cluster1 <- names(markers1)[i]
    print(cluster1)
    cluster1_markers <- markers1[[i]]
    
    for(j in 1:length(markers2)){
      cluster2 <- names(markers2)[j]
      print(cluster2)
      cluster2_markers <- markers2[[j]]
      
      similarity <- similarityHypergeometric(cluster1_markers, cluster2_markers, background)
      similarityJaccard(markers1, markers2, background, TRUE)
      print(similarity)
      scores[n] <- similarity
      clusters1[n] <- cluster1
      clusters2[n] <- cluster2
      n <- n + 1
    }
  }
  scores <- data.frame("markers1" = clusters1, "markers2" = clusters2, "Hypergeometric_pval" = scores)
  scores$negLogPval <- -1 * log10(scores$Hypergeometric_pval)
  names(scores)[names(scores) == "markers1"] <- name_markers1
  names(scores)[names(scores) == "markers2"] <- name_markers2
  name_score <- "negLogPval"
  print(head(scores))
  g <- ggplot(scores, aes_string(name_markers2, name_score, fill = name_markers1))
  g <- g + geom_col()
  g <- g + geom_hline(yintercept = -1 * log10(0.05))
  g <- g + facet_grid(reformulate(name_markers1, "."), labeller = label_both)
  g + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none")
}

# get all possible human genes from our GATA6 experiment:
load("joshua_workspace_data/GATA6-mmB-Clustering 7-18.RData")
GATA6_all_genes <- row.names(mmB_D5)
GATA6_markers <- mmB_D5.markers
n_GATA6_cluster <- length(unique(GATA6_markers$cluster))
GATA6_markers_list <- vector(mode = "list", length = n_GATA6_cluster)
names(GATA6_markers_list) <- unique(GATA6_markers$cluster)
for(i in 1:length(GATA6_markers_list)){
  cluster <- names(GATA6_markers_list)[i]
  GATA6_markers_list[[i]] <- GATA6_markers$gene[GATA6_markers$cluster == cluster]
}
rm(mmB_D5)
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(s.genes)
rm(g2m.genes)

# get all possible mouse genes from the mouse endoderm explorer
nowotschin <- readRDS("nowotschin_data/nowotschin_E6.5_filtered_seurat_object.RDS")
nowotschin_all_genes <- row.names(nowotschin)
rm(nowotschin)
nowotschin_E6.5_markers <- readRDS("nowotschin_data/nowotschin_E6.5_filtered_seurat_markers.RDS")
n_nowotschin_E6.5_clusters <- length(unique(nowotschin_E6.5_markers$cluster))
nowotschin_markers_list <- vector(mode = "list", length = n_nowotschin_E6.5_clusters)
names(nowotschin_markers_list) <- unique(nowotschin_E6.5_markers$cluster)
for(i in 1:length(nowotschin_markers_list)){
  cluster <- names(nowotschin_markers_list)[i]
  markers <- nowotschin_E6.5_markers$gene[nowotschin_E6.5_markers$cluster == cluster]
  nowotschin_markers_list[[i]] <- convertMouseGeneList(markers)
}


# create background set as human genes
nowotschin_human_orthologs <- convertMouseGeneList(nowotschin_all_genes)
background <- intersect(GATA6_all_genes, nowotschin_human_orthologs)

g <- compareMarkersBetweenDatasets(GATA6_markers_list, nowotschin_markers_list, background, "GATA6_cluster", "Nowotschin_E6.5_celltype")
g + labs(x="Nowotschin E6.5 Cell Type", y="-log10(Hyperg Pval)")

compareMarkersWithinDataset(GATA6_markers_list)

compareMarkersWithinDataset(nowotschin_markers_list)


# read in markers from nowotschin E4.5
# read in markers from nowotschin E6.5
# read in markers from sala