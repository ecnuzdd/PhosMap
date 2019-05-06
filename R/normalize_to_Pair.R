#' For data with pairs, normalize them to the sample with flag eaqul to -1.
#'
#' @param data_frame a data frame as input.
#' @param experiment_design_file a data frame for design of experiment.
#'
#' @return A data frame after normalization.
#' @export
normalize_to_Pair <- function(data_frame, experiment_design_file){
  data_frame_colnames <- colnames(data_frame)
  ID <- as.vector(data_frame[,1])
  Value_raw <- data_frame[,-1]
  Value_raw_colnames <- colnames(Value_raw)
  groups_labels <- names(table(experiment_design_file$Group))
  groups_labels_len <- length(groups_labels)
  correction_df <- NULL
  correction_df_colnames <- NULL
  for(i in seq_len(groups_labels_len)){
    group_label <- groups_labels[i]
    case_index <- which(experiment_design_file$Group == group_label & experiment_design_file$Pair == 1)
    case_v <- as.vector(unlist(Value_raw[,case_index]))
    control_index <- which(experiment_design_file$Group == group_label & experiment_design_file$Pair == (-1))
    control_v <- as.vector(unlist(Value_raw[,control_index]))
    correction_v <- (case_v + 1) / (control_v + 1)
    correction_df <- cbind(correction_df, correction_v)
    correction_df_colname <- paste(group_label, Value_raw_colnames[case_index], Value_raw_colnames[control_index], sep = '_')
    correction_df_colnames <- c(correction_df_colnames, correction_df_colname)
  }

  data_frame_normalization <- data.frame(ID, correction_df)
  colnames(data_frame_normalization) <- c(data_frame_colnames[1], correction_df_colnames)
  return(data_frame_normalization)
}
