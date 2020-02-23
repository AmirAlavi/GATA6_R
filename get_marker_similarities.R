library(ggplot2)

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

# get all possible human genes from our GATA6 experiment:
load("joshua_workspace_data/GATA6-mmB-Clustering 7-18.RData")
GATA6_all_genes <- row.names(mmB_D5)
GATA6_markers <- mmB_D5.markers
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

# create background set as human genes
nowotschin_human_orthologs <- convertMouseGeneList(nowotschin_all_genes)
background <- intersect(GATA6_all_genes, nowotschin_human_orthologs)




nowotschin_markers_as_hgnc <- vector(mode = "list", length = n_nowotschin_E6.5_clusters)
names(nowotschin_markers_as_hgnc) <- unique(nowotschin_E6.5_markers$cluster)
for(i in 1:length(nowotschin_markers_as_hgnc)){
  cluster <- names(nowotschin_markers_as_hgnc)[i]
  print(cluster)
  markers <- nowotschin_E6.5_markers$gene[nowotschin_E6.5_markers$cluster == cluster]
  nowotschin_markers_as_hgnc[[i]] <- convertMouseGeneList(markers)
}

n_GATA6_cluster <- length(unique(GATA6_markers$cluster))
n_scores <- 2 * n_GATA6_cluster^2
jaccard_sims <- numeric(n_scores)
gata6_cluster1 <- numeric(n_scores)
gata6_cluster2 <- numeric(n_scores)
part <- character(n_scores)
# analyze the similarity of the within-dataset cluster markers
i <- 1
for(cluster1 in unique(GATA6_markers$cluster)) {
  markers1 <- GATA6_markers$gene[GATA6_markers$cluster == cluster1]
  for(cluster2 in unique(GATA6_markers$cluster)) {
    markers2 <- GATA6_markers$gene[GATA6_markers$cluster == cluster2]
    jaccard <- similarityJaccard(markers1, markers2, verbose = TRUE)
    gata6_cluster1[i] <- cluster1
    gata6_cluster2[i] <- cluster2
    jaccard_sims[i] <- jaccard
    part[i] <- "similarity"
    
    gata6_cluster1[i+1] <- cluster1
    gata6_cluster2[i+1] <- cluster2
    jaccard_sims[i+1] <- 1 - jaccard
    part[i+1] <- "remainder"
    
    i <- i + 2
  }
}

scores <- data.frame("GATA6_clusterA" = gata6_cluster1, "GATA6_clusterB" = gata6_cluster2, "Jaccard_Similarity" = jaccard_sims, "Partition" = part)
g <- ggplot(scores, aes(x=1, y=Jaccard_Similarity, fill = Partition))
g <- g + geom_col(color = 'black')
g <- g + coord_polar("y", start=0)
g <- g + facet_grid(rows = vars(GATA6_clusterB), cols = vars(GATA6_clusterA))
g <- g + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
g
# g + theme_minimal()
# g + theme_void()
#g <- g + labs(x="Nowotschin E6.5 Cell Type", y="-log10(Hyperg Pval)")
#g + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none")


compareMarkerLists <- function(markers1, markers2, background, name_markers1, name_markers2, species1 = c("human", "mouse"), species2 = c("human", "mouse")) {
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
      
      similarity <- similarityHypergeometric(markers1, markers2, background)
      scoress[n] <- similarity
      clusters1[n] <- cluster1
      clusters2[n] <- cluster2
      n <- n + 1
    }
  }
  scores <- data.frame(name_markers1 = clusters1, name_markers2 = clusters2, "Hypergeometric_pval" = scores)
  g <- ggplot(scores, aes(name_markers2, -1 * log10(Hypergeometric_pval), fill = name_markers1))
  g <- g + geom_col()
  g <- g + geom_hline(yintercept = -1 * log10(0.05))
  g <- g + facet_grid(cols = vars(name_markers1), labeller = label_both)
  g + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none")
}
g <- compareMarkerLists()
g + labs(x="Nowotschin E6.5 Cell Type", y="-log10(Hyperg Pval)")

n_nowotschin_E6.5_clusters <- length(unique(nowotschin_E6.5_markers$cluster))
n_scores = n_GATA6_cluster * n_nowotschin_E6.5_clusters
jaccard_sims <- numeric(n_scores)
hyperg_sims <- numeric(n_scores)
gata6_clusters <- numeric(n_scores)
nowotschin_celltypes <- character(n_scores)
i <- 1
for(GATA6_cluster in unique(GATA6_markers$cluster)){
  print(GATA6_cluster)
  GATA6_cluster_markers <- GATA6_markers$gene[GATA6_markers$cluster == GATA6_cluster]
  
  for(j in 1:length(nowotschin_markers_as_hgnc)){
    nowots_clust <- names(nowotschin_markers_as_hgnc)[j]
    nowots_markers <- nowotschin_markers_as_hgnc[[j]]
    jaccard_sim <- similarityJaccard(GATA6_cluster_markers, nowots_markers, background, TRUE)
    print(jaccard_sim)
    hyperg_sim <- similarityHypergeometric(GATA6_cluster_markers, nowots_markers, background, TRUE)
    print(hyperg_sim)
    gata6_clusters[i] <- GATA6_cluster
    nowotschin_celltypes[i] <- nowots_clust
    jaccard_sims[i] <- jaccard_sim
    hyperg_sims[i] <- hyperg_sim
    i <- i + 1
  }
}

scores <- data.frame("GATA6_cluster" = gata6_clusters, "Nowotschin_E6.5_celltype" = nowotschin_celltypes, "Jaccard_Similarity" = jaccard_sims, "Hypergeometric_pval" = hyperg_sims)

neg_log10_trans <- scales::trans_new("neg_log10", transform = function(x) -1 * log10(x), inverse = function(x) 10^(-x), domain = c(.Machine$double.xmin, Inf))

g <- ggplot(scores, aes(Nowotschin_E6.5_celltype, -1 * log10(Hypergeometric_pval), fill = GATA6_cluster))
g <- g + geom_col()
g <- g + geom_hline(yintercept = -1 * log10(0.05))
# t <- t + scale_y_continuous(trans = neg_log10_trans)
g <- g + facet_grid(cols = vars(GATA6_cluster), labeller = label_both)
g <- g + labs(x="Nowotschin E6.5 Cell Type", y="-log10(Hyperg Pval)")
g + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position = "none")



# read in markers from nowotschin E4.5
# read in markers from nowotschin E6.5
# read in markers from sala