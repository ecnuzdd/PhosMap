#' For data without pairs but with control, normalize them to the control.
#'
#' @param data_frame a data frame as input.
#' @param experiment_design_file a data frame for design of experiment.
#' @param control_label a string for a control.
#'
#' @return A data frame after normalization.
#' @export

normalize_nopair_ctrl_by_col <- function(data_frame, experiment_design_file, control_label){
  data_frame_colnames <- colnames(data_frame)
  ID <- as.vector(data_frame[,1])
  Value_raw <- data_frame[,-1]

  control_index <- which(grepl(control_label, experiment_design_file$Group))
  control_Value <- Value_raw[,control_index]
  if(length(control_index) == 1){
    control_Value <- as.data.frame(control_Value)
  }
  control_Value_row_mean <- apply(control_Value, 1, mean)
  Value_new <- (Value_raw + 1)/(control_Value_row_mean + 1)
  data_frame_normalization <- data.frame(ID, Value_new)
  colnames(data_frame_normalization) <- data_frame_colnames
  return(data_frame_normalization)

}
