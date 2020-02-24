library(Seurat)
library(plyr)
source("utils.R")
# alignment

# Clusters:
#   0, 1, 2 exVE
# 4: Epi
# 5: emVE
# 3: maybe mesoderm, but probably not in E6.5 (more likely to be germs cell like)

load("joshua_workspace_data/GATA6-mmB-Clustering 7-18.RData")
rm(mmB_D5.data)
rm(mmB_D5.markers)
rm(g2m.genes)
rm(s.genes)
GATA6_all_genes <- row.names(mmB_D5)




nowotschin <- readRDS("nowotschin_data/nowotschin_E6.5_filtered_seurat_object_processed.RDS")
nowotschin_all_genes <- row.names(nowotschin)

HGNC2MGI_map <- getHuman2MouseGeneMapping(GATA6_all_genes)

HGNC2MGI_map <- HGNC2MGI_map[!duplicated(HGNC2MGI_map$HGNC.symbol), ]
HGNC2MGI_map <- HGNC2MGI_map[!duplicated(HGNC2MGI_map$MGI.symbol), ]

common_genes <- HGNC2MGI_map[is.element(HGNC2MGI_map$MGI.symbol, nowotschin_all_genes), ]

# dt[is.element(dt$fct, vc),]
# 
# 
# nowotschin_human_orthologs <- convertMouseGeneList(nowotschin_all_genes)
# common_genes <- intersect(GATA6_all_genes, nowotschin_human_orthologs)
# common_human2mouse_map <- getHuman2MouseGeneMapping(common_genes)
# common_human2mouse_map <- common_human2mouse_map[!duplicated(common_human2mouse_map$HGNC.symbol),]
# common_human2mouse_map <- common_human2mouse_map[!duplicated(common_human2mouse_map$MGI.symbol), ]



human <- subset(x = mmB_D5, features = common_genes$HGNC.symbol)
mouse <- subset(x = nowotschin, features = common_genes$MGI.symbol)


integrated <- FindIntegrationAnchors(object.list = c(human, mouse), dims = 1:30)
