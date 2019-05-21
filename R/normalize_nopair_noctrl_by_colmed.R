#' For data without pairs and control, normalize them to the median.
#'
#' @param data_frame a data frame as input.
#'
#' @return A data frame after normalization.
#' @export
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/normalize_nopair_noctrl_by_colmed.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "normalize_nopair_noctrl_by_colmed.RData")
#' load("normalize_nopair_noctrl_by_colmed.RData")
#'
#' phospho_data_normalize_by_column <- normalize_nopair_noctrl_by_colmed(
#'   phospho_data_normalized
#' )
#' head(phospho_data_normalize_by_column)
#' }


normalize_nopair_noctrl_by_colmed <- function(data_frame){
  requireNamespace('stats')
  data_frame_colnames <- colnames(data_frame)
  ID <- as.vector(data_frame[,1])
  Value_raw <- data_frame[,-1]
  Value_median <- Value_raw
  Value_median_row <- nrow(Value_median)
  for(i in seq_len(Value_median_row)){
    x <- as.vector(unlist(Value_raw[i,]))
    x_median <- stats::median(x)
    Value_median[i,] <- (x+1)/(x_median+1)
  }
  data_frame_normalization <- data.frame(ID, Value_median)
  colnames(data_frame_normalization) <- data_frame_colnames
  return(data_frame_normalization)
}
