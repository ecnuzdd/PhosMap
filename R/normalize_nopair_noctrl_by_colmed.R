#' For data without pairs and control, normalize them to the median.
#'
#' @param data_frame a data frame as input.
#'
#' @return A data frame after normalization.
#' @export
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
