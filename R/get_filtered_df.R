#' Get data frame filtered based on the Mascot and reference files.
#'
#'
#' @param mascotfileName a string for mascot names as input.
#' @param refFileName a string for reference file names.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @return A filtered data frame
#'
#' @export

get_filtered_df <- function(mascotfileName, refFileName){
  requireNamespace('utils')
  fileData <- utils::read.table(mascotfileName, header = TRUE, sep = '\t')
  fileData_seq <- as.vector(fileData$pep_seq)

  refFileData <- utils::read.delim(refFileName)
  refFileData_seq <- as.vector(refFileData$Sequence)

  fileData_seq_upper <- toupper(fileData_seq)
  refFileData_seq_upper <- toupper(refFileData_seq)
  refFileData_seq_upper <- refFileData_seq_upper[which(!duplicated(refFileData_seq_upper))]
  refFileData_seq_upper_len <- length(refFileData_seq_upper)
  kept_index <- NULL
  for(i in seq_len(refFileData_seq_upper_len)){
    tmp_seq_upper <- refFileData_seq_upper[i]
    index <- which(fileData_seq_upper==tmp_seq_upper)
    kept_index <- c(kept_index, index)
  }

  fileData_subset <- fileData[kept_index,]
  fileData_subset_query_num <- fileData_subset$pep_query_num
  fileData_subset_query_num_unique <- unique(fileData_subset_query_num)
  fileData_subset_query_num_unique_len <- length(fileData_subset_query_num_unique)
  index1 <- NULL
  for(i in seq_len(fileData_subset_query_num_unique_len)){
    tmp_query_num <- fileData_subset_query_num_unique[i]
    index <- which(fileData_subset_query_num==tmp_query_num)
    tmp_query_rank <- fileData_subset$query_rank[index]
    tmp_query_rank_max <- min(tmp_query_rank)
    tmp_query_rank_max_index <- which(tmp_query_rank_max==tmp_query_rank)
    tmp_index <- index[tmp_query_rank_max_index]
    index1 <- c(index1, tmp_index)
  }
  reuslt_df <- fileData_subset[index1,]

  reuslt_df_seq <- as.vector(reuslt_df$pep_seq)
  reuslt_df_seq_unique <- unique(reuslt_df_seq)
  reuslt_df_seq_unique_len <- length(reuslt_df_seq_unique)
  index2 <- NULL
  for(i in seq_len(reuslt_df_seq_unique_len)){
    tmp_seq_unique <- reuslt_df_seq_unique[i]
    index <- which(reuslt_df_seq==tmp_seq_unique)
    tmp_pep_score <- reuslt_df$pep_score[index]
    pep_score_max_index <- which.max(tmp_pep_score)
    kept_index <- index[pep_score_max_index]
    index2 <- c(index2, kept_index)
  }

  reuslt_df1 <- reuslt_df[index2,]
  pep_var_mods <- as.vector(reuslt_df1$pep_var_mod)
  pep_var_mod_confs <- as.vector(reuslt_df1$pep_var_mod_conf)
  delete_index <- which(pep_var_mods!='' & pep_var_mod_confs=='')

  reuslt_df2 <- reuslt_df1[-delete_index,]

  return(reuslt_df2)
}
