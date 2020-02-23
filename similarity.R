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

