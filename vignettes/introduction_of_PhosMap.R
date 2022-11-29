## ----load---------------------------------------------------------------------
# load PhosMap
library('PhosMap')
# load intermediate results from https://github.com/ecnuzdd/PhosMap_datasets
# manully download from https://github.com/ecnuzdd/PhosMap_datasets/blob/master/data/BRAFi.RData
load("BRAFi.RData")

# background_df A data frame for motif enrichment analysis as background
# combined_df_with_mapped_gene_symbol Get a data frame mapped GI number to Gene Symbol
# data_frame_normalization_with_control_no_pair A data frame containning phosphoproteomics data normalized by proteomics data.
# foreground_df A data frame for motif enrichment analysis as foreground.
# fuzzy_input_df A data frame for time course analysis as input.
# group A factor for experiment group information.
# merge_df_with_phospho_peptides A merged phosphoproteomics data frame based on peptides files (unique ID).
# motif_group_m_ratio_df_mat A matrix for motif profile.
# phospho_data_filtering_STY_and_normalization A phosphoproteomics data frame after normalization and filtering.
# profiling_data_normalized A proteomics data frame after normalization and filtering.
# summary_df_of_unique_proteins_with_sites A data frame that phosphorylation sites had been mapping to protein sequence and eliminated redundancy.


## ----extract_psites_score, eval = FALSE---------------------------------------
#  BASE_DIR <- getwd() # working directory
#  BASE_DIR <- normalizePath(BASE_DIR)
#  phosphorylation_exp_design_info_file_path <- normalizePath(file.path(BASE_DIR, 'phosphorylation_exp_design_info.txt'))
#  phosphorylation_peptide_dir <- normalizePath(file.path(BASE_DIR, 'phosphorylation_peptide_txt'))
#  
#  if(FALSE){
#    # if you have xml files from mascot results, you can run the cmd to parser them to text files.
#    mascot_xml_dir <- normalizePath(file.path(BASE_DIR, 'mascot_xml'))
#    mascot_txt_dir <- normalizePath(file.path(BASE_DIR, 'mascot_txt'))
#    extract_psites_score(phosphorylation_exp_design_info_file_path, mascot_xml_dir, mascot_txt_dir)
#  
#    # Based on above-mentioned text files from Mascot results,
#    # the following cmd can generate CSV files of phosphorylation sites with confidence score.
#    psites_score_dir <- normalizePath(file.path(BASE_DIR, 'psites_score_txt'))
#    generate_psites_score_file(mascot_txt_dir, phosphorylation_peptide_dir, psites_score_dir)
#  }
#  # Merge phosphoproteomics data based on peptides files (unique ID).
#  # If qc = TRUE, considering confidence score of phosphorylation sites.
#  # A merged phosphoproteomics data frame based on peptides files (unique ID).
#  merge_df_with_phospho_peptides <- pre_process_filter_psites(
#    phosphorylation_peptide_dir,
#    psites_score_dir,
#    phosphorylation_exp_design_info_file_path,
#    qc = TRUE,
#    min_score = 20,
#    min_FDR = 0.01
#  )

## -----------------------------------------------------------------------------
# Get a data frame mapped GI number to Gene Symbol.
# system.time({
#   combinated_df_with_mapped_gene_symbol = get_combined_data_frame(
#     merge_df_with_phospho_peptides, species = 'human', id_type = 'RefSeq_Protein_GI'
#   )
# })
head(combined_df_with_mapped_gene_symbol)

## -----------------------------------------------------------------------------
# Assign psites to protein sequence.
# Unique ID: protein_gi + phosphorylation site in protein sequence.
# system.time({
#   summary_df_of_unique_proteins_with_sites = get_summary_with_unique_sites(
#     combinated_df_with_mapped_gene_symbol, species = 'human', fasta_type = 'refseq'
#   )
# })
head(summary_df_of_unique_proteins_with_sites)

## ----eval = FALSE-------------------------------------------------------------
#  # Imputation with the next order of magnitude of the minimum except for zero.
#  # Filtering data only including phosphorylation site.
#  phospho_data_filtering_STY_and_normalization_list <- get_normalized_data_of_psites(
#    summary_df_of_unique_proteins_with_sites,
#    phosphorylation_exp_design_info_file_path,
#    topN = NA, mod_types = c('S', 'T', 'Y')
#  )
#  phospho_data_filtering_STY <- phospho_data_filtering_STY_and_normalization_list$ptypes_area_df_with_id
#  phospho_data_filtering_STY_and_normalization <- phospho_data_filtering_STY_and_normalization_list$ptypes_fot5_df_with_id
#  head(phospho_data_filtering_STY_and_normalization)

## ----eval = FALSE-------------------------------------------------------------
#  # Based on phospho_data_filtering_STY_and_normalization
#  ID <- paste(phospho_data_filtering_STY_and_normalization$GeneSymbol,
#             phospho_data_filtering_STY_and_normalization$AA_in_protein,
#             sep = '_')
#  Value <- phospho_data_filtering_STY_and_normalization[,-seq(1,6)]
#  phospho_data <- data.frame(ID, Value)
#  phospho_data_rownames <- paste(phospho_data_filtering_STY_and_normalization$GI,
#                                phospho_data_filtering_STY_and_normalization$GeneSymbol,
#                                phospho_data_filtering_STY_and_normalization$AA_in_protein,
#                                sep = '_')
#  rownames(phospho_data) <- phospho_data_rownames
#  # Further normalize phosphoproteomics data based on proteomics data
#  # The configurations of function see help document.
#  data_frame_normalization_with_control_no_pair <- normalize_phos_data_to_profiling(
#    phospho_data, profiling_data_normalized,
#    phosphorylation_exp_design_info_file_path,
#    profiling_exp_design_info_file_path,
#    control_label = '0',
#    pair_flag = FALSE
#  )
#  head(data_frame_normalization_with_control_no_pair)

## ----visualization_with_simple_tsne, fig.cap = "clustering example: t-SNE", fig.wide = TRUE----
# Clustering: t-SNE or PCA
expr_data_frame <- data_frame_normalization_with_control_no_pair
# t-SNE using all experiments
visualization_with_simple_tsne(expr_data_frame, group)

## ----visualization_with_simple_pca, fig.cap = "clustering example: PCA", fig.wide = TRUE----
expr_ID <- as.vector(expr_data_frame[,1])
expr_Valule <- expr_data_frame[,-1]
expr_Valule_mean <- NULL
expr_Valule_row <- nrow(expr_Valule)
for(i in 1:expr_Valule_row){
  x <- as.vector(unlist(expr_Valule[i,]))
  x_m <- tapply(x, group, mean)
  expr_Valule_mean <- rbind(expr_Valule_mean, x_m)
}
group_levels = levels(group)
colnames(expr_Valule_mean) <- group_levels
expr_df <- data.frame(expr_ID, expr_Valule_mean)
# PCA using mean value in group for comparison with original literature
visualization_with_simple_pca(expr_df)

## -----------------------------------------------------------------------------
# Differently expressed Proteins/Genes analysis
# t2 vs t0
expr_data_frame <- data_frame_normalization_with_control_no_pair[,1:17] # phosphoproteomics data normalized by proteomics data
# select phosphorylation sites with greater variation
expr_data_frame_var <- apply(expr_data_frame, 1, function(x){
  var(x[-1])
})
index_of_kept <- which(expr_data_frame_var>1)
expr_data_frame <- expr_data_frame[index_of_kept,]

# group information (t0 vs t2)
deps_group_levels <- c('t0', 't2')
deps_group <- factor(as.vector(group)[1:16], levels = deps_group_levels)

## ----analysis_deps_limma, fig.cap = "Differential expression analysis: limma", fig.wide = TRUE----
# (1) limma
limma_results_df <- analysis_deps_limma(expr_data_frame, deps_group, deps_group_levels, log2_label = FALSE, adjust_method = 'none')
limma_results_df$ID <- apply(limma_results_df, 1, function(x){
  x = strsplit(x, '_')[[1]]
  paste(x[2], x[3], sep = '_')
})
visualization_deps_with_scatter(limma_results_df, minFC = 2, minPvalue = 0.05, main = 'Differentially expressed proteins  \n with limma',
                                show_text = FALSE, min_up_text = 70, min_down_text = 70)

## ----eval=FALSE---------------------------------------------------------------
#  # (2) SAM
#  sam_results_list <- analysis_deps_sam(expr_data_frame, deps_group, log2_label = FALSE, nperms = 100, rand = NULL, minFDR = 0.05, samr_plot = TRUE)
#  sam_results <- rbind(sam_results_list$genes_up_df, sam_results_list$genes_down_df)

## ----eval=FALSE---------------------------------------------------------------
#  # (3) annova
#  anova_result_df <- analysis_deps_anova(expr_data_frame, deps_group, log2_label = FALSE, return_padjust = TRUE, adjust_method = 'BH')
#  visualization_deps_with_scatter(anova_result_df, minFC = 2, minPvalue = 0.05, main = 'Differentially expressed proteins \n with anova',
#                                  show_text = FALSE, min_up_text = 15, min_down_text = 15)

## ----visualization_fuzzycluster, fig.cap = "Time course analysis example", fig.width = 6, fig.height = 6, fig.align = 'center'----
group_levels <- levels(group)
# fuzzy c-means clustering
set.seed(1000)
fuzzy_clustObj <- visualization_fuzzycluster(
	fuzzy_input_df, group, group_levels, 
	k_cluster=9, iteration = 100, 
	mfrow = c(3,3), min_mem = 0.1,
	plot = TRUE
)
# clusters information
clusterS_info <- fuzzy_clustObj$cluster
clusterS_names <- names(clusterS_info)
clusters_df <- data.frame(clusterS_names, clusterS_info)
# write.csv(clusters_df, 'clusters_df.csv', row.names = TRUE)

## -----------------------------------------------------------------------------
# For early and late response
# early -> clusterS_info==1
# late -> clusterS_info==2
cluster_flag <- 'early'
cluster_symbol <- clusterS_names[clusterS_info==1]
expr_data_frame <- data_frame_normalization_with_control_no_pair
index_of_cluster <- match(cluster_symbol, expr_data_frame$ID)
cluster_df <- expr_data_frame[index_of_cluster,]

