#' For data without pairs but with control, normalize them to the control.
#'
#' @param data_frame a data frame as input.
#' @param experiment_design_file a data frame for design of experiment.
#' @param control_label a string for a control.
#'
#' @return A data frame after normalization.
#' @export
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url <- "https://github.com/ecnuzdd/PhosMap_datasets/function_demo_data/normalize_nopair_ctrl_by_col.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'Rdata')
#' writeBin(load_data, "normalize_nopair_ctrl_by_col.RData")
#' load("normalize_nopair_ctrl_by_col.RData")
#'
#' phospho_data_normalize_by_column <- normalize_nopair_ctrl_by_col(
#'   phospho_data_normalized,
#'   phosphorylation_experiment_design_file,
#'   control_label
#' )
#' head(phospho_data_normalize_by_column)
#' }



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
