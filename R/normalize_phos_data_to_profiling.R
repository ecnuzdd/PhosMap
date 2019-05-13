#' Normalize phospho-data to profiling
#'
#' @param phospho_data_normalized A data frame of phospho-data after normalization
#' @param profiling_data_normalized A data frame of profiling after normalization
#' @param phosphorylation_exp_design_info_file_path A file path about phosphorylation experiment design, it has 2 kinds of file configuration as follows:
#'        1. experiment_design_noPair.txt must contain columns of Experiment_Code, Group.
#'        2. experiment_design_Pair.txt must contain columns of Experiment_Code, Group, and Pair. (Pair: 1 -> case, -1 -> control)
#'
#' @param profiling_exp_design_info_file_path A file path about profiling experiment design, it has 2 kinds of file configuration as same as phosphorylation_exp_design_info_file_path.
#'
#' @param control_label A string represents label of control group. The default is NA which shows no control group.
#' @param pair_flag A boolean value represents whether experiments have pairs. The default is FALSE which shows no pairs.
#'
#' @return A data frame which comes from results that phospho-data is normalizated base on the abundance of proteins in the profiling experiments.
#' @export
#'
#' @examples
#' demo_data_url1 <- url('https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/function_demo_data/phospho_data_topX.RData')
#' demo_data_url2 <- url('https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/function_demo_data/profiling_data_normalized.RData')
#' demo_data_url3 <- url('https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/function_demo_data/phosphorylation_exp_design_info.txt')
#' demo_data_url4 <- url('https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/function_demo_data/profiling_exp_design_info.txt')
#'
#' load(demo_data_url1)
#' load(demo_data_url2)
#' phosphorylation_exp_design_info_file_path <- demo_data_url3
#' profiling_exp_design_info_file_path <- demo_data_url4
#'
#' data_frame_normalization_with_control_no_pair <- normalize_phos_data_to_profiling(
#'   phospho_data_topX, profiling_data_normalized,
#'   phosphorylation_exp_design_info_file_path,
#'   profiling_exp_design_info_file_path,
#'   control_label = '0',
#'   pair_flag = FALSE)
#'
#' head(data_frame_normalization_with_control_no_pair)
#'

normalize_phos_data_to_profiling <- function (phospho_data_normalized, profiling_data_normalized,
                                                 phosphorylation_exp_design_info_file_path, profiling_exp_design_info_file_path,
                                                 control_label = NA, pair_flag = FALSE) {

  requireNamespace('utils')
  # phosphorylation_exp_design_info_file_path
  phosphorylation_experiment_design_file <- utils::read.table(phosphorylation_exp_design_info_file_path, header = TRUE,
                                                       sep = '\t', stringsAsFactors = NA)
  phosphorylation_groups_labels <- names(table(phosphorylation_experiment_design_file$Group))
  phosphorylation_groups <- factor(phosphorylation_experiment_design_file$Group,
                                  levels = phosphorylation_groups_labels)


  # profiling_exp_design_info_file_path
  profiling_experiment_design_file <- utils::read.table(profiling_exp_design_info_file_path,
                                                header = TRUE, sep = '\t', stringsAsFactors = NA)
  profiling_groups_labels <- names(table(profiling_experiment_design_file$Group))
  profiling_groups <- factor(profiling_experiment_design_file$Group, levels = profiling_groups_labels)


  if(!pair_flag){
    if(is.na(control_label)){ # no control, nomalization with median
      cat('\n', 'No pair and no control, nomalization with median')
      # normalize phospho-data
      phospho_data_normalize_by_column <- normalize_nopair_noctrl_by_colmed(phospho_data_normalized)
      # normalize profiling
      profiling_data_normalizedby_column <- normalize_nopair_noctrl_by_colmed(profiling_data_normalized)
    }else{ # having control, nomalization with median
      cat('\n', 'No pair and having control, nomalization with median of control')
      phospho_data_normalize_by_column <- normalize_nopair_ctrl_by_col(phospho_data_normalized, phosphorylation_experiment_design_file, control_label)
      profiling_data_normalizedby_column <- normalize_nopair_ctrl_by_col(profiling_data_normalized, profiling_experiment_design_file, control_label)
    }
  }else{
    cat('\n', 'Having pair, case is nomalized to control')
    phospho_data_normalize_by_column <- normalize_to_Pair(phospho_data_normalized, phosphorylation_experiment_design_file)
    profiling_data_normalizedby_column <- normalize_to_Pair(profiling_data_normalized, profiling_experiment_design_file)


  }

  phospho_ID <- as.vector(phospho_data_normalize_by_column[,1])
  phospho_ID_Symbol <- apply(data.frame(phospho_ID), 1, function(x){
    x <- strsplit(x, split = '_')[[1]]
    x[1] # gene symbol
  })
  phospho_Value <- (phospho_data_normalize_by_column[,-1])
  phospho_Value_vs_profiling <- phospho_Value
  phospho_Value_row <- nrow(phospho_Value)
  phospho_Value_col <- ncol(phospho_Value)

  profiling_ID <- as.vector(profiling_data_normalizedby_column[,1])
  profiling_Value <- profiling_data_normalizedby_column[,-1]

  # phosphorylation_groups_labels
  # profiling_groups_labels
  if(pair_flag){
    for(i in seq_len(phospho_Value_row)){
      x_phospho_Symbol <- phospho_ID_Symbol[i]
      x_phospho <- as.vector(unlist(phospho_Value[i,]))
      j <- which(profiling_ID == x_phospho_Symbol)
      if(length(j) > 0){
        x_correction <- NULL
        x_profiling <- as.vector(unlist(profiling_Value[j,]))
        for(k in seq_len(phospho_Value_col)){
          phos_k_group_label <- phosphorylation_groups_labels[k]
          prof_k_vector <- x_profiling[which(phos_k_group_label==profiling_groups_labels)]
          prof_k_mean <- mean(prof_k_vector)
          x_k_correction <- x_phospho[k]/prof_k_mean
          x_correction <- c(x_correction, x_k_correction)
        }
      }else{
        x_correction <- x_phospho
      }
      phospho_Value_vs_profiling[i,] <- x_correction
    }
  }else{
    for(i in seq_len(phospho_Value_row)){
      x_phospho_Symbol <- phospho_ID_Symbol[i]
      x_phospho <- as.vector(unlist(phospho_Value[i,]))
      j <- which(profiling_ID == x_phospho_Symbol)
      if(length(j) > 0){
        x_correction <- NULL
        x_profiling <- as.vector(unlist(profiling_Value[j,]))
        for(k in seq_len(phospho_Value_col)){
          phos_k_group_label <- phosphorylation_experiment_design_file$Group[k]
          prof_k_vector <- x_profiling[which(phos_k_group_label==profiling_experiment_design_file$Group)]
          prof_k_mean <- mean(prof_k_vector)
          x_k_correction <- x_phospho[k]/prof_k_mean
          x_correction <- c(x_correction, x_k_correction)
        }
      }else{
        x_correction <- x_phospho
      }
      phospho_Value_vs_profiling[i,] <- x_correction
    }
  }

  df_phospho_Value_vs_profiling <- data.frame(phospho_ID, phospho_Value_vs_profiling)
  colnames(df_phospho_Value_vs_profiling) <- colnames(phospho_data_normalize_by_column)
  cat('\n', 'Normalization over. (phosphorylation vs profiling)')
  return(df_phospho_Value_vs_profiling)

}
