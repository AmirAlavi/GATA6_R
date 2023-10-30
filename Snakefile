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
        "results/figures/data_integration_and_alignment/mmB_D5_onto_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/mmB_D5_onto_Tyser_transfer_label_distribution.png",
        "results/figures/data_integration_and_alignment/D4_merged_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D4_merged_Tyser_transfer_label_distribution.png",
        "results/figures/data_integration_and_alignment/D4_SB431521_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D4_SB431521_Tyser_transfer_label_distribution.png",
        "results/figures/data_integration_and_alignment/D5_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D5_Tyser_transfer_label_distribution.png",
        "results/figures/data_integration_and_alignment/D7_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D7_Tyser_transfer_label_distribution.png"

###############################################################################
# Data downloading and marker finding
###############################################################################
rule download_our_data:
    output:
        "remote_data/Ours/GATA6-mmB-Clustering 7-18.RData"
    run:
        gdown.download(GDRIVE_URLS["Our_data"], output[0])

rule ensure_local_data:
    output:
        "data/iDiscoid_Seurat_Objects/mmBmK12hr_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmK_36hr_03302022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmK_3hr_03302022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0D1_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0D5_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0Ri_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD1_030722_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD4_combined_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD5_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD7_processed.RData"
    run:
        print("\n********************\nERROR: Missing iDiscoid Seurat Objects, download them to:\n\tdata/iDiscoid_Seurat_Objects\n********************\n")

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

rule create_xiang_new_annotations_marker_lists_by_day:
    input:
        "results/R_objects/Xiang_new_annotations_processed_seurat_object.RDS"
    output:
        "results/R_objects/Xiang_new_annotations_markers_by_day.RDS",
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/find_xiang_new_annotation_markers_by_day.R"

rule assemble_our_revision_data_markers_by_day:
    input:
        "data/iDiscoid_Seurat_Objects/mmBmK12hr_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmK_36hr_03302022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmK_3hr_03302022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0D1_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0D5_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD0Ri_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD1_030722_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD4_combined_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD5_033022_processed.RData",
        "data/iDiscoid_Seurat_Objects/mmBmKD7_processed.RData"
    output:
        "results/R_objects/iDiscoid_markers_by_day.RDS",
        "results/R_objects/iDiscoid_gene_sets_by_day.RDS"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/assemble_our_markers_by_day.R"
        

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
        #metadata="remote_data/Tyser_et_al/umap.rds"
        metadata="remote_data/Tyser_et_al/tyser_annot_umap.rds"
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

rule revision_comparisons:
    input:
        "results/R_objects/Xiang_new_annotations_processed_seurat_object.RDS",
        "results/R_objects/Xiang_new_annotations_markers_by_day.RDS",
        "results/R_objects/iDiscoid_markers_by_day.RDS",
        "results/R_objects/iDiscoid_gene_sets_by_day.RDS",
        "results/R_objects/Tyser_markers.RDS",
        "results/R_objects/Tyser_seurat_object_processed.RDS",
    output:
        "results/figures/hypergeometric_comparisons/revision/success.txt"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_xiang_and_tyser_revision.R"

rule aug22_revision_comparison:
    input:
        Tyser_seurat_obj = "results/R_objects/Tyser_seurat_object_processed.RDS",
        Tyser_markers = "results/R_objects/Tyser_markers.RDS",
        Ma_counts = "remote_data/Ma_et_al/Counts Matrix - GSE130114_MF1453.csv",
        Ma_markers = "remote_data/Ma_et_al/aax7890-Ma-SM-Table-S6.csv",
        Nowotschin_E6pt5_seurat_obj = "results/R_objects/Nowotschin_E6.5_seurat_object_processed.RDS",
        Nowotschin_E4pt5_markers = "data/Nowotschin_et_al/E4.5_markers.csv",
        Nowotschin_E5pt5_markers = "data/Nowotschin_et_al/E5.5_markers.csv",
        Nowotschin_E6pt5_markers = "results/R_objects/Nowotschin_E6.5_markers.RDS",
        Nowotschin_E7pt5_markers = "data/Nowotschin_et_al/E7.5_markers.csv",
        Sala_markers = "data/Sala_et_al/Sala_markers.csv",
        Sala_all_genes = "data/Sala_et_al/genes.tsv",
        Xiang_seurat_obj = "results/R_objects/Xiang_new_annotations_processed_seurat_object.RDS",
        Xiang_markers = "results/R_objects/Xiang_new_annotations_markers_by_day.RDS",
        mmBmK_12hr = "data/iDiscoid_Seurat_Objects/mmBmK12hr_033022_processed.RData",
        mmBmK_36hr = "data/iDiscoid_Seurat_Objects/mmBmK_36hr_03302022_processed.RData",
        mmBmK_3hr = "data/iDiscoid_Seurat_Objects/mmBmK_3hr_03302022_processed.RData",
        mmBmK_D0D1 = "data/iDiscoid_Seurat_Objects/mmBmKD0D1_033022_processed.RData",
        mmBmK_D0D5 = "data/iDiscoid_Seurat_Objects/mmBmKD0D5_033022_processed.RData",
        mmBmK_D0Ri = "data/iDiscoid_Seurat_Objects/mmBmKD0Ri_033022_processed.RData",
        mmBmK_D1 = "data/iDiscoid_Seurat_Objects/mmBmKD1_030722_processed.RData",
        mmBmK_D2 = "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        mmBmK_D3 = "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData",
        mmBmK_D4_combined = "data/iDiscoid_Seurat_Objects/mmBmKD4_combined_processed_WTsubclustered.RData",
        mmBmK_D5 = "data/iDiscoid_Seurat_Objects/mmBmKD5_033022_processed.RData",
        mmBmK_D7 = "data/iDiscoid_Seurat_Objects/mmBmKD7_processed.RData",
        D4_merged = "data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData",
        D4_SB431542_Aug22Revisions = "data/Aug22_revision_seurat_objects/D4_SB431521_Aug22Revision.RData",
        D5_20220330_Aug22Revisions = "data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData",
        D7_20220307_Aug22Revisions = "data/Aug22_revision_seurat_objects/D7_Aug22Revision.RData",
        FeLO_monoculture_cytokine_minus = "data/Aug22_revision_seurat_objects/FeLO-monoculture-cytokine-minus_03152022.RData",
        FeLO_monoculture_cytokine_plus = "data/Aug22_revision_seurat_objects/FeLO-monoculture-cytokine-plus_03152022.RData",
        D8_H1iDiscoid_Aug22Revision = "data/Aug22_revision_seurat_objects/H1iDiscoidD8_033022_Aug22Revision.RData",
        yolk_sac = "data/Aug22_revision_seurat_objects/Yolk_sac_subclustered.RData"
    output:
        outfile = "results/figures/hypergeometric_comparisons/success_aug22rev.txt"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_Aug22Revision.R"
        
rule fetalLiver_yolksac_comparison:
    input:
        "data/fetalLiver_yolksac_comparison_objects/coculture_FeLO_03152022.RData",
        "data/fetalLiver_yolksac_comparison_objects/FeLO-monoculture-cytokine-minus_03152022.RData",
        "data/fetalLiver_yolksac_comparison_objects/FeLO-monoculture-cytokine-plus_03152022.RData",                       "data/fetalLiver_yolksac_comparison_objects/Yolk_sac_subclustered.RData",
        "data/fetalLiver_yolksac_comparison_objects/YolkSacFetalLiver-split days with markers.RData",                     "data/fetalLiver_yolksac_comparison_objects/AdultPlusFetalLiverYolkSac-split days with markers.RData",
        "data/fetalLiver_yolksac_comparison_objects/ZhaiCynomolgusEmbryoProcessed-split days with markers.RData",
        "data/fetalLiver_yolksac_comparison_objects/DesLO_D17_4-24-23.rds",
        "data/fetalLiver_yolksac_comparison_objects/FeLO_D17_JH-04-20-23.rds",
        "data/fetalLiver_yolksac_comparison_objects/YSFetalLiverFullMetadata-markers.rds",
        "data/fetalLiver_yolksac_comparison_objects/DesLO_D17_5-03-23-markers.rds",
        "data/fetalLiver_yolksac_comparison_objects/FeLO_D17_5-03-23-markers.rds"
    output:
        "results/figures/hypergeometric_comparisons/success_fetalLiver_yolksac.txt"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/hyperg_similarity_fetalLiver_yolkSac_comparison.R"
        
rule plot_hypergeom_comparison:
    input:
        "results/figures/hypergeometric_comparisons/success_fetalLiver_yolksac.txt",
        "results/figures/hypergeometric_comparisons/success_aug22rev.txt"
    output:
        "results/figures/hypergeometric_comparisons/success_comp_plot.txt"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/plot_hypergeom.R"

###############################################################################
# YSEndo VS DE(P) analysis
###############################################################################
rule ysendo_dep_analysis:
    input:
        D4_merged = "data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData",
        mmBmK_D2 = "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        D5_20220330_Aug22Revisions = "data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData",
        mmBmK_D3 = "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData"
    output:
        outfile = "results/figures/comparison_idiscoid_to_YSEndo_and_DE/success_ysendo_dep_plot.txt"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/comparison_idiscoid_to_YSEndo_and_DE.R"

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
rule integrate_and_align_mmB_D5_tyser:
    input:
        "results/R_objects/Ours_subclusters_seurat_object.RDS",
        # "results/R_objects/Tyser_markers.RDS",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/tyser_annot_umap.rds"
    output:
        "results/figures/data_integration_and_alignment/mmB_D5_onto_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/mmB_D5_onto_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_mmB_D5-tyser.R"
        
rule integrate_and_align_D4_merged_tyser:
    input:
        "data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/tyser_annot_umap.rds"
    output:
        "results/figures/data_integration_and_alignment/D4_merged_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D4_merged_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_D4_merged-tyser.R"

rule integrate_and_align_D4_SB431521_Aug22Revision_tyser:
    input:
        "data/Aug22_revision_seurat_objects/D4_SB431521_Aug22Revision.RData",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/tyser_annot_umap.rds"
    output:
        "results/figures/data_integration_and_alignment/D4_SB431521_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D4_SB431521_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_D4_SB431521-tyser.R"

rule integrate_and_align_D5_tyser:
    input:
        "data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/tyser_annot_umap.rds"
    output:
        "results/figures/data_integration_and_alignment/D5_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D5_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_D5-tyser.R"

rule integrate_and_align_D7_tyser:
    input:
        "data/Aug22_revision_seurat_objects/D7_Aug22Revision.RData",
        "remote_data/Tyser_et_al/express_vals.rds",
        "remote_data/Tyser_et_al/tyser_annot_umap.rds"
    output:
        "results/figures/data_integration_and_alignment/D7_Tyser_UMAP.png",
        "results/figures/data_integration_and_alignment/D7_Tyser_transfer_label_distribution.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/alignment_D7-tyser.R"

###############################################################################
# Generate expression change plots separately for WT and iGATA6 clusters in D0~D5
###############################################################################
rule plot_expression_change_D0_to_D5:
    input:
        mmBmK_D0Ri = "data/iDiscoid_Seurat_Objects/mmBmKD0Ri_033022_processed.RData",
        mmBmK_12hr = "data/iDiscoid_Seurat_Objects/mmBmK12hr_033022_processed.RData",
        mmBmK_D1 = "data/iDiscoid_Seurat_Objects/mmBmKD1_030722_processed.RData",
        mmBmK_36hr = "data/iDiscoid_Seurat_Objects/mmBmK_36hr_03302022_processed.RData",
        mmBmK_D2 = "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        mmBmK_D3 = "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData",
        D4_merged = "data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData",
        D5_20220330_Aug22Revisions = "data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData"
    output:
        WT = "results/figures/expression_change/WT_expression_change.png",
        GATA6 = "results/figures/expression_change/iGATA6_expression_change.png"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/expression_change_D0_to_D5.R"

###############################################################################
# Generate inputs for scdiff2
###############################################################################
rule generate_input_for_scdiff2:
    input:
        mmBmK_D0Ri = "data/iDiscoid_Seurat_Objects/mmBmKD0Ri_033022_processed.RData",
        mmBmK_12hr = "data/iDiscoid_Seurat_Objects/mmBmK12hr_033022_processed.RData",
        mmBmK_D1 = "data/iDiscoid_Seurat_Objects/mmBmKD1_030722_processed.RData",
        mmBmK_36hr = "data/iDiscoid_Seurat_Objects/mmBmK_36hr_03302022_processed.RData",
        mmBmK_D2 = "data/iDiscoid_Seurat_Objects/mmBmKD2_033022_processed.RData",
        mmBmK_D3 = "data/iDiscoid_Seurat_Objects/mmBmKD3_03302022_processed_WTsubclustered.RData",
        D4_merged = "data/Aug22_revision_seurat_objects/D4_Aug22Revision_combined.RData",
        D5_20220330_Aug22Revisions = "data/Aug22_revision_seurat_objects/D5_Aug22Revision.RData",
        cluster_label = "data/iDiscoid_cluster_assignments/D0_to_D5_Cluster_Assignments.csv"
    output:
        scdiff2_out_WT = "results/scdiff2/expression_WT",
        scdiff2_out_GATA6 = "results/scdiff2/expression_GATA6"
    conda:
        "envs/r_seurat_env.yml"
    script:
        "scripts/generate_input_for_scdiff2.R"