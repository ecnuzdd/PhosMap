#' To normalize data and filter data only including phosphorylation sites.
#'
#' @param data_frame A data frame containing IDs and quantification values merged from multi-experiments as input.
#' @param experiment_code_file_path A file path of storing experiment codes as input. The experiment codes are required to keep pace with column names of Value.
#' @param topN, A numeric value, selecting N p-sites with high intensity rank for normalization, the default is NA.
#' @param mod_types, A vector for modification residues, the default is c('S', 'T', 'Y') for phosphorylation modifications.
#'
#' @return A list including data frame after filtering or normalization (x 1e5).
#' @export
#'
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url1 <- "https://github.com/ecnuzdd/PhosMap_datasets/function_demo_data/get_normalized_data_of_psites.RData"
#' ftp_url2 <- "https://github.com/ecnuzdd/PhosMap_datasets/function_demo_data/phosphorylation_exp_design_info.txt"

#' load_data1 <- load_data_with_ftp(ftp_url1, 'Rdata')
#' writeBin(load_data1, "get_normalized_data_of_psites.RData")
#' load("get_normalized_data_of_psites.RData")
#'
#' load_data2 <- load_data_with_ftp(ftp_url2, 'downloadtxt')
#' writeBin(load_data2, "phosphorylation_exp_design_info.txt")
#' phosphorylation_exp_design_info_file_path <- "./phosphorylation_exp_design_info.txt"
#'
#' phospho_data_filtering_STY_and_normalization_list <- get_normalized_data_of_psites(
#'   summary_df_of_unique_proteins_with_sites,
#'   phosphorylation_exp_design_info_file_path,
#'   topN = NA, mod_types = c('S', 'T', 'Y')
#' )
#' head(phospho_data_filtering_STY_and_normalization_list)
#'
#' }


get_normalized_data_of_psites <- function(data_frame, experiment_code_file_path, topN = NA, mod_types = c('S', 'T', 'Y')){
  requireNamespace('utils')
  cat('\n The 7th step: Normalize data and filter data only including phosphorylation site.')
  experiment_code <- utils::read.table(experiment_code_file_path, header = TRUE, sep = '\t', stringsAsFactors = NA)
  experiment_code <- as.vector(unlist(experiment_code$Experiment_Code))
  data_frame_colnames <- colnames(data_frame)

  cat('\n The 7th step is running.')
  summary_df_ID_Info <- data_frame[, seq_len(6)]
  summary_df_ID_Info$AA_in_protein <- toupper(summary_df_ID_Info$AA_in_protein)
  summary_df_Value <- data_frame[, -(seq_len(6))]

  cat('\n Filtering data only including S/T/Y modifications.')
  ptypes <- mod_types
  index_of_AA_in_protein <- apply(data.frame(summary_df_ID_Info$AA_in_protein), 1, function(x){
    if(grepl('S', x) | grepl('T', x) | grepl('Y', x)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  index_of_ptypes <- which(index_of_AA_in_protein)
  if(length(index_of_ptypes)>0){
    ptypes_id_df <- summary_df_ID_Info[index_of_ptypes,]
    ptypes_value <- summary_df_Value[index_of_ptypes,]
  }else{
    message('No data with modifications taking place on ', paste(mod_types, collapse = '|'))
    stop('')
  }

  cat('\n Normalization based on the total sum.')
  Value_FOT5 <- ptypes_value
  Value_FOT5_col <- ncol(Value_FOT5)
  if(is.na(topN)){
    for(i in seq_len(Value_FOT5_col)){
      x <- as.vector(unlist(ptypes_value[,i]))
      Value_FOT5[,i] <- x/sum(x)*1e5
    }
  }else{
    for(i in seq_len(Value_FOT5_col)){
      x <- as.vector(unlist(ptypes_value[,i]))
      x_order <- order(x, decreasing = TRUE)
      x_order_top <- x_order[seq_len(topN)]
      x[-x_order_top] <- 0
      Value_FOT5[,i] <- x/sum(x)*1e5
    }
  }
  ptypes_value_FOT5 <- as.matrix(Value_FOT5)
  cat('\n Imputation with the next order of magnitude of the minimum except for zero.')
  index_of_zero <- which(ptypes_value_FOT5==0)
  min_value_of_non_zero <- min(ptypes_value_FOT5[-index_of_zero])
  ptypes_value_FOT5[index_of_zero] <- min_value_of_non_zero*0.1

  ptypes_df_list <- list(
    ptypes_area_df_with_id = data.frame(ptypes_id_df, ptypes_value),
    ptypes_fot5_df_with_id = data.frame(ptypes_id_df, ptypes_value_FOT5)
  )

  cat('\n The 7th step is over ^_^.')
  return(ptypes_df_list)
}
