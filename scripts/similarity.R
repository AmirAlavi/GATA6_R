library(ggplot2)

source("scripts/utils.R")

#' Compute the similarity between two sets using a Hypergeometric test for overrepresentation.
#' Equivalent to a one-tailed Fisher's Exact Test.
#'
#' @param a A set of items (e.g. list of genes)
#' @param b Another set of items (e.g. list of genes)
#' @param background Background set of all possible items that could have been in \code{a} or \code{b}. I.e. the "universe" of items
#' @param verbose Whether or not to print intermediate values
#'
#' @return The pvalue result of the hypergeometric test
#' @export
#'
#' @examples pval <- similarityHypergeometric(gene_list1, gene_list2, only_consider_these_genes_list)
similarityHypergeometric <- function(a, b, background = NULL, verbose = FALSE){
  a <- unique(a)
  b <- unique(b)
  if (!is.null(background)) {
    background <- unique(background)
    a <- intersect(a, background)
    b <- intersect(b, background)
    
    n_urn <- length(background)
  } else {
    n_urn <- union(a, b)
  }
  
  n_white_balls <- length(a)
  n_black_balls <- n_urn - n_white_balls
  n_draws <- length(b)
  n_drawn_successes <- length(intersect(a, b))
  if (verbose) {
    print(paste("\t(similarityHypergeometric) n_urn:", n_urn))
    print(paste("\t(similarityHypergeometric) n_white_balls:", n_white_balls))
    print(paste("\t(similarityHypergeometric) n_draws:", n_draws))
    print(paste("\t(similarityHypergeometric) n_draws_white:", n_drawn_successes))
  }
  
  phyper(n_drawn_successes-1, n_white_balls, n_black_balls, n_draws, lower.tail = FALSE)
}

#' Compute the Jaccard similarity between two sets 
#'
#' @param a A set of items (e.g. list of genes)
#' @param b Another set of items (e.g. list of genes)
#' @param background Background set of all possible items that could have been in \code{a} or \code{b}. I.e. the "universe" of items
#' @param verbose Whether or not to print intermediate values
#'
#' @return The jaccard similarity between \code{a} and \code{b}
#' @export
#'
#' @examples similarity <- similarityJaccard(gene_list1, gene_list2, only_consider_these_genes_list)
similarityJaccard <- function(a, b, background = NULL, verbose = FALSE){
  a <- unique(a)
  b <- unique(b)
  if (!is.null(background)) {
    background <- unique(background)
    a <- intersect(a, background)
    b <- intersect(b, background)
  }
  
  intersect_ab <- intersect(a, b)
  union_ab <- union(a, b)
  if (verbose) {
    print(paste("\t(similarityJaccard) intersect:", length(intersect_ab)))
    print(paste("\t(similarityJaccard) union:", length(union_ab)))
  }
  
  length(intersect_ab) / length(union_ab)
}

#' Plot the Jaccard similarity between marker sets within a dataset
#'
#' @param markers A named list where each element is a list of markers for a cluster/phenotype, and the corresponding name is the cluster/phenotype name
#'
#' @return A ggplot2 object containing a graph, represented as a matrix-like figure where each element is a radial plot showing the Jaccard similarity between the column marker set and row marker set
#' @export
#'
#' @examples
compareMarkersWithinDataset <- function(markers) {
  n_scores <- 2 * length(markers) ^ 2
  jaccard_sims <- numeric(n_scores)
  clusters1 <- numeric(n_scores)
  clusters2 <- numeric(n_scores)
  part <- character(n_scores)
  
  n <- 1
  for (i in 1:length(markers)) {
    cluster1 <- names(markers)[i]
    cluster1_markers <- markers[[i]]
    for (j in 1:length(markers)) {
      cluster2 <- names(markers)[j]
      cluster2_markers <- markers[[j]]
      
      jaccard <- similarityJaccard(cluster1_markers, cluster2_markers)
      
      clusters1[n] <- cluster1
      clusters2[n] <- cluster2
      jaccard_sims[n] <- jaccard
      part[n] <- "similarity"
      
      clusters1[n + 1] <- cluster1
      clusters2[n + 1] <- cluster2
      jaccard_sims[n + 1] <- 1 - jaccard
      part[n + 1] <- "remainder"
      
      n <- n + 2
    }
  }
  
  scores <- data.frame("A" = clusters1, "B" = clusters2, "Jaccard_Similarity" = jaccard_sims, "Partition" = part)
  g <- ggplot(scores, aes(x = 1, y = Jaccard_Similarity, fill = Partition))
  g <- g + geom_col(color = 'black')
  g <- g + coord_polar("y", start = 0)
  g <- g + facet_grid(rows = vars(B), cols = vars(A))
  g + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
}

#' Plot the Hypergeometric similarity pvalues between the marker sets of two datasets
#'
#' @param markers1 A named list where each element is a list of markers for a cluster/phenotype, and the corresponding name is the cluster/phenotype name
#' @param markers2 A named list where each element is a list of markers for a cluster/phenotype, and the corresponding name is the cluster/phenotype name
#' @param background Background set to use in Hypergeometric test
#' @param name_markers1 String name of the dataset for \code{markers1}
#' @param name_markers2 String name of the dataset for \code{markers2}
#'
#' @return A ggplot2 object containing a graph, showing the similarities of each cluster in \code{markers1} with each cluster in \code{markers2} as a bar chart, higher bars -> lower hypergeometric pvals -> more similarity.
#' @export
#'
#' @examples
compareMarkersBetweenDatasets <- function(markers1, markers2, background, name_markers1, name_markers2) {
  n_scores = length(markers1) * length(markers2)
  scores <- numeric(n_scores)
  clusters1 <- character(n_scores)
  clusters2 <- character(n_scores)
  
  n <- 1
  for (i in 1:length(markers1)) {
    cluster1 <- names(markers1)[i]
    cluster1_markers <- markers1[[i]]
    
    for (j in 1:length(markers2)) {
      cluster2 <- names(markers2)[j]
      cluster2_markers <- markers2[[j]]
      
      similarity <- similarityHypergeometric(cluster1_markers, cluster2_markers, background)
      
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
  #g <- ggplot(scores, aes_string(name_markers2, name_score, fill = name_markers1))
  #g <- g + geom_col()
  #g <- g + geom_hline(yintercept = -1 * log10(0.05))
  #pal <- getColorPalette(names(markers1))
  #g <- g + scale_fill_manual(values = pal)
  #g <- g + facet_grid(reformulate(name_markers1, "."), labeller = label_both)
  #g + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
  return(scores)
}
