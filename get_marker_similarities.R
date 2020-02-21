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

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}