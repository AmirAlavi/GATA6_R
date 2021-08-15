"""
Requires the 'snakemake' conda environment (defined in environment.yml):
$ conda activate snakemake

Then, to run all analysis and generate all figures, run:
$ snakemake --cores 1 --use-conda

To generate only a particular result, just specify it, e.g.:
$ snakemake --cores 1 --use-conda results/figures/hypergeometric_comparisons/Nowotschin_E6.5.png
"""
import gdown

# Large data stored in google drive. The files below are all from this folder:
# https://drive.google.com/drive/folders/1JQqgd82EKzcIYOIGTD9YkI5IVIX73fZz?usp=sharing
GDRIVE_URLS = {
    "Our_data": "https://drive.google.com/uc?id=1Bk-_i4Pcwi3ijNpcm89JuIhBOE15qG-F",
    "Nowotschin_counts": "https://drive.google.com/uc?id=1MQBDdtmyMas7YikQXjhFWrsvZony1HZd",
    "Nowotschin_metadata": "https://drive.google.com/uc?id=1C6oOEDjRwTcz9Txs-UEhFGmBYihr4SXw",
    "Xiang_counts": "https://drive.google.com/uc?id=1MtF7qyWfI_GgBEqS7747fPD4j08SXhqt",
    "Xiang_metadata": "https://drive.google.com/uc?id=12e3tV6S84WX1Um7_YPBTyADZ5xoWOyAO",
    "Tyser_counts": "https://drive.google.com/uc?id=1PiPLHNG2eUIMI31OkgVZE33DBhC43Nmc",
    "Tyser_metadata": "https://drive.google.com/uc?id=1N0Ap1KFElx-IDhqIwBTY8kkh2QI4gRIl",
    "Ma_counts": "https://drive.google.com/uc?id=1xbiuxtumjD9t_OG-Dd7toxRt6g5e5Ub5",
    "Ma_markers": "https://drive.google.com/uc?id=1KGuC9uMRgZQAvtXDnll7HzSnPZaS-oi2", 
    "Sala_gene_list": "https://drive.google.com/uc?id=1fIq8MzTsFln0l87S-8GXZEuqjvjNZc54",
}


rule all:
    input:
        "results/figures/hypergeometric_comparisons/Nowotschin_E4.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E5.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_disjoint.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E7.5.png",
        "results/figures/hypergeometric_comparisons/Ours_intra_marker_similarity.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_intra_marker_similarity.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_disjoint_intra_marker_similarity.png",
        "results/figures/hypergeometric_comparisons/Tyser.png",
        "results/figures/hypergeometric_comparisons/Xiang.png",
        "results/figures/hypergeometric_comparisons/Tyser_and_Xiang.png",
        "results/figures/hypergeometric_comparisons/Xiang_new_annotations.png",
        "results/figures/hypergeometric_comparisons/Sala.png",
        "results/figures/hypergeometric_comparisons/Ma.png",
        "results/figures/marker_module_scoring/Xiang_new_annotations_TSNE.png",
        "results/figures/marker_module_scoring/Xiang_new_annotations_UMAP.png",
        "results/figures/data_integration_and_alignment/Ours_onto_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/Ours_onto_Tyser_transfer_label_distribution.png"

###############################################################################
# Data downloading and marker finding
###############################################################################
rule download_our_data:
    output:
        "remote_data/Ours/GATA6-mmB-Clustering 7-18.RData"
    run:
        gdown.download(GDRIVE_URLS["Our_data"], output[0])

rule subcluster_our_data:
    input:
        "remote_data/Ours/GATA6-mmB-Clustering 7-18.RData"
    output:
        "results/R_objects/Ours_subclusters_seurat_object.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/subcluster.R"

rule download_nowotschin_data:
    output:
        counts="remote_data/Nowotschin_et_al/sc_endoderm_all_cells_counts.csv",
        metadata="remote_data/Nowotschin_et_al/sc_endoderm_all_cells_metadata.csv"
    run:
        gdown.download(GDRIVE_URLS["Nowotschin_counts"], output.counts)
        gdown.download(GDRIVE_URLS["Nowotschin_metadata"], output.metadata)

rule create_nowotschin_seurat_object:
    input:
        "remote_data/Nowotschin_et_al/sc_endoderm_all_cells_counts.csv",
        "remote_data/Nowotschin_et_al/sc_endoderm_all_cells_metadata.csv"
    output:
        "results/R_objects/Nowotschin_E6.5_seurat_object.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/create_nowotschin_seurat_object.R"

rule create_nowotschin_marker_lists:
    input:
        "results/R_objects/Nowotschin_E6.5_seurat_object.RDS"
    output:
        "results/R_objects/Nowotschin_E6.5_markers.RDS",
        "results/R_objects/Nowotschin_E6.5_seurat_object_processed.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/find_nowotschin_E6.5_markers.R"

rule download_xiang_data:
    output:
        counts="remote_data/Xiang_et_al/counts.csv",
        metadata="remote_data/Xiang_et_al/meta.csv"
    run:
        gdown.download(GDRIVE_URLS["Xiang_counts"], output.counts)
        gdown.download(GDRIVE_URLS["Xiang_metadata"], output.metadata)
    
rule create_xiang_seurat_object:
    input:
        "remote_data/Xiang_et_al/counts.csv",
        "remote_data/Xiang_et_al/meta.csv"
    output:
        "results/R_objects/Xiang_seurat_object.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/create_xiang_seurat_object_counts.R"

rule create_xiang_new_annotations_marker_lists:
    input:
        "results/R_objects/Xiang_seurat_object.RDS",
        "data/Xiang_et_al_counts/Xiang_merged_label_table_analysisColumn.csv"
    output:
        "results/R_objects/Xiang_new_annotations_markers.RDS",
        "results/R_objects/Xiang_new_annotations_processed_seurat_object.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/find_xiang_new_annotation_markers.R"

rule download_tyser_data:
    output:
        counts="remote_data/Tyser_et_al/express_vals.rds",
        metadata="remote_data/Tyser_et_al/umap.rds"
    run:
        gdown.download(GDRIVE_URLS["Tyser_counts"], output.counts)
        gdown.download(GDRIVE_URLS["Tyser_metadata"], output.metadata)

rule create_tyser_marker_lists:
    input:
        counts="remote_data/Tyser_et_al/express_vals.rds",
        metadata="remote_data/Tyser_et_al/umap.rds"
    output:
        "results/R_objects/Tyser_markers.RDS",
        "results/R_objects/Tyser_seurat_object_processed.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/find_tyser_markers.R"
        
rule download_ma_data:
    output:
        counts="remote_data/Ma_et_al/Counts Matrix - GSE130114_MF1453.csv",
        markers="remote_data/Ma_et_al/aax7890-Ma-SM-Table-S6.csv"
    run:
        gdown.download(GDRIVE_URLS["Ma_counts"], output.counts)
        gdown.download(GDRIVE_URLS["Ma_markers"], output.markers)


rule download_sala_gene_list:
    output:
        "data/Sala_et_al/genes.tsv.gz"
    run:
        gdown.download(GDRIVE_URLS["Sala_gene_list"], output[0])

rule extract_sala_gene_list:
    input:
        "data/Sala_et_al/genes.tsv.gz"
    output:
        "data/Sala_et_al/genes.tsv"
    shell:
        "gunzip data/Sala_et_al/genes.tsv.gz"

###############################################################################
# Hypergeometric similarity tests
###############################################################################
rule compare_to_nowotschin:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        "results/R_objects/Nowotschin_E6.5_seurat_object_processed.RDS",
        "results/R_objects/Nowotschin_E6.5_markers.RDS"
    output:
        "results/figures/hypergeometric_comparisons/Nowotschin_E4.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E5.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_disjoint.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E7.5.png",
        "results/figures/hypergeometric_comparisons/Ours_intra_marker_similarity.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_intra_marker_similarity.png",
        "results/figures/hypergeometric_comparisons/Nowotschin_E6.5_disjoint_intra_marker_similarity.png",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_nowotschin.R"

rule compare_to_tyser_and_xiang:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        "results/R_objects/Tyser_markers.RDS",
        "results/R_objects/Tyser_seurat_object_processed.RDS",
        "results/R_objects/Xiang_seurat_object.RDS",
        "results/R_objects/Xiang_new_annotations_markers.RDS",
        "results/R_objects/Xiang_new_annotations_processed_seurat_object.RDS"
    output:
        "results/figures/hypergeometric_comparisons/Tyser.png",
        "results/figures/hypergeometric_comparisons/Xiang.png",
        "results/figures/hypergeometric_comparisons/Tyser_and_Xiang.png",
        "results/figures/hypergeometric_comparisons/Xiang_new_annotations.png",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_xiang_and_tyser.R"

rule compare_to_sala:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        "data/Sala_et_al/genes.tsv"
    output:
        "results/figures/hypergeometric_comparisons/Sala.png",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_sala.R"

rule compare_to_ma:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        "remote_data/Ma_et_al/Counts Matrix - GSE130114_MF1453.csv",
        "remote_data/Ma_et_al/aax7890-Ma-SM-Table-S6.csv"
    output:
        "results/figures/hypergeometric_comparisons/Ma.png",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_ma.R"

###############################################################################
# Marker "module score" figures
###############################################################################
rule create_module_score_figures:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        "results/R_objects/Xiang_new_annotations_markers.RDS",
    output:
        "results/figures/marker_module_scoring/Xiang_new_annotations_TSNE.png",
        "results/figures/marker_module_scoring/Xiang_new_annotations_UMAP.png",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/module_score_figures.R"

###############################################################################
# Data integration and alignment figures
###############################################################################
rule integrate_and_align_tyser:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        # "results/R_objects/Tyser_markers.RDS",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/umap.rds"
    output:
        "results/figures/data_integration_and_alignment/Ours_onto_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/Ours_onto_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_tyser.R"