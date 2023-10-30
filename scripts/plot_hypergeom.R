library(ggplot2)
library(Seurat)

source("scripts/utils.R")

input_dir <- "results/figures/hypergeometric_comparisons/figure_data"
output_dir <- "results/figures/hypergeometric_comparisons/figures"
score_files <- list.files(path = input_dir, pattern = "*.csv")
out_files <- gsub(".csv$","",score_files)
all_pvals <- c()
all_scores <- list()

# Read hyperg results
for(i in 1:length(score_files)){
    scores <- read.csv(paste(input_dir, score_files[i],sep = "/"))
    all_pvals <- c(all_pvals,scores$Hypergeometric_pval)
    all_scores[[length(all_scores)+1]]<-scores
}

# Adjust p-values
all_pvals_adj <- p.adjust(all_pvals, method = "BH")

# Assign adjusted p-values
start = 1
for(i in 1:length(all_scores)){
    end = start + nrow(all_scores[[i]])-1
    all_scores[[i]]$Hypergeometric_pval <- all_pvals_adj[start:end]
    start = end + 1
}

dir.create('results/figures/hypergeometric_comparisons/figures')
output_dir <- 'results/figures/hypergeometric_comparisons/figures'

# Generate all plots with adjusted p-values
for(i in 1:length(all_scores)){
    scores <- all_scores[[i]]
    name_markers1 <- colnames(scores)[1]
    name_markers2 <- colnames(scores)[2]
    scores[[name_markers1]] <- as.character(scores[[name_markers1]])
    scores[[name_markers2]] <- as.character(scores[[name_markers2]])
    g <- ggplot(scores, aes_string(name_markers2, "negLogPval", fill = name_markers1))
    g <- g + geom_col()
    g <- g + geom_hline(yintercept = -1 * log10(0.05))
    pal <- getColorPalette(unique(scores[[name_markers1]]))
    g <- g + scale_fill_manual(values = pal)
    g <- g + facet_grid(reformulate(name_markers1, "."), labeller = label_both)
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
    g <- g + labs(x = paste(name_markers2,"Cell Type"), y = "-log10(Hyperg Pval)")
    
    # Save plot and plot data
    outname <- paste0(output_dir,
                      "/",
                       out_files[i],
                       ".png"
                      )
    ggsave(outname, plot=g, width = 10*length(unique(scores[[name_markers1]]))/12, height = 5, units = "in") # For better visualization, automatic adjust the width so that labels won't overlap with each other.
}

file.create(snakemake@output[[1]])