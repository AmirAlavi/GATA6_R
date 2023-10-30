# Basic functions to convert between mouse and human gene names
# Adapted from https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/

#' Get a mapping of HGNC symbols to MGI symbols (human to mouse)
#'
#' @param x A list or vector of HGNC query gene symbols to map
#'
#' @return A 2 column dataframe with the mappings
#' @export
#'
#' @examples h2m_df <- getMouse2HumanGeneMapping(hgnc_genes)
getMouse2HumanGeneMapping <- function(x) {
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = x ,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = T
  )
}

#' Get a mapping of MGI symbols to HGNC symbols (mouse to human)
#'
#' @param x A list or vector of MGI query gene symbols to map
#'
#' @return A 2 column dataframe with the mappings
#' @export
#'
#' @examples m2h_df <- getHuman2MouseGeneMapping(mgi_genes)
getHuman2MouseGeneMapping <- function(x) {
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x ,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = T
  )
}

#' Convert MGI symbols to HGNC symbols (mouse to human)
#'
#' @param x A list or vector of MGI query gene symbols to convert
#' 
#' @param verbose Whether to print preview of first few hits
#'
#' @return A list with the converted genes, not necessarily the same length or order as query
#' @export
#'
#' @examples human_orthos <- convertMouseGeneList(mgi_genes)
convertMouseGeneList <- function(x, verbose = FALSE) {
  genesV2 <- getMouse2HumanGeneMapping(x)
  humanx <- unique(genesV2[, 2])
  
  if (verbose) {
    # Print the first 6 genes found to the screen
    print(head(humanx))
  }
  return(humanx)
}

#' Convert HGNC symbols to MGI symbols (human to mouse)
#'
#' @param x A list or vector of HGNC query gene symbols to convert
#' 
#' @param verbose Whether to print preview of first few hits
#'
#' @return A list with the converted genes, not necessarily the same length or order as query
#' @export
#'
#' @examples mouse_orthos <- convertHumanGeneList(hgnc_genes)
convertHumanGeneList <- function(x, verbose = FALSE) {
  genesV2 <- getHuman2MouseGeneMapping(x)
  mousex <- unique(genesV2[, 2])
  
  if (verbose) {
    # Print the first 6 genes found to the screen
    print(head(mousex))
  }
  return(mousex)
}


# Consistent Color Palette
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

getColorPalette <- function(idents, shift = 9) {
  pal <- DiscretePalette(26, "alphabet2")
  pal <- shifter(pal, 9)
  n.unused.colors <- length(pal) - length(idents)
  idents <- c(idents, rep("none", n.unused.colors))
  names(pal) <- idents
  return(pal)
}