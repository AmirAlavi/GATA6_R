# list_of_markers , each element is a vector of markers, each element also has a name
# list_of_markers2

# for each list in lm1
# for each list in lm2
# similarity (lm1, lm2)

# human GATA6 data
# mouse endoderm data

# get all the genes from the mouse data
# convert them to a list of orthologous human genes
# subset those to the intersect with the human GATA6 genes
# https://support.bioconductor.org/p/96955/

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

# get all possible mouse genes from the mouse endoderm explorer
nowotschin <- readRDS("nowotschin_data/nowotschin_E6.5_filtered_seurat_object.RDS")
nowotschin_all_genes <- row.names(nowotschin)
rm(nowotschin)
nowotschin_E6.5_markers <- readRDS("nowotschin_data/nowotschin_E6.5_filtered_seurat_markers.RDS")

# create background set as human genes
nowotschin_human_orthologs <- convertMouseGeneList(nowotschin_all_genes)
background <- intersect(GATA6_all_genes, nowotschin_human_orthologs)



for(GATA6_cluster in unique(GATA6_markers$cluster)){
  print(GATA6_cluster)
  GATA6_cluster_markers <- GATA6_markers$gene[GATA6_markers$cluster == GATA6_cluster]
  for(nowotschin_E6.5_cluster in unique(nowotschin_E6.5_markers$cluster)){
    print(nowotschin_E6.5_cluster)
    nowotschin_E6.5_cluster_markers <- nowotschin_E6.5_markers$gene[nowotschin_E6.5_markers$cluster == nowotschin_E6.5_cluster]
    nowotschin_E6.5_cluster_markers <- convertMouseGeneList(nowotschin_E6.5_cluster_markers)
    jaccard_sim <- similarityJaccard(GATA6_cluster_markers, nowotschin_E6.5_cluster_markers, background, TRUE)
    hyperg_sim <- similarityHypergeometric(GATA6_cluster_markers, nowotschin_E6.5_cluster_markers, background, TRUE)
    print(jaccard_sim)
    print(hyperg_sim)
  }
}

# read in markers from nowotschin E4.5
# read in markers from nowotschin E6.5
# read in markers from sala