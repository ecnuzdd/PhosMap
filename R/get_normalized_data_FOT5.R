#' Normailization on basis of sum
#'
#' @param data_frame A data frame containing IDs and values merged from multi-experiments as input.
#' @param experiment_code_file_path A file path of storing experiment codes as input. The experiment codes are required to keep pace with column names of Values.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame after normalization
#' @export
#'
#' @examples
#' \dontrun{
#' data_frame_normalization = get_normalized_data_FOT5(
#'   data_frame,
#'   experiment_code_file_path
#' )
#' }
get_normalized_data_FOT5 <- function(
  data_frame,
  experiment_code_file_path
){
  requireNamespace('utils')
  cat('\n The 7th step: Normalize data and filter data only including phosphorylation site.')
  experiment_code <- utils::read.table(experiment_code_file_path, header = TRUE, sep = '\t', stringsAsFactors = NA)
  experiment_code <- as.vector(unlist(experiment_code$Experiment_Code))
  data_frame_colnames <- colnames(data_frame)
  ID <- as.vector(data_frame[,1])
  Value_raw <- data_frame[,-1]
  Value_FOT5 <- Value_raw
  Value_FOT5_col <- ncol(Value_FOT5)
  for(i in seq_len(Value_FOT5_col)){
    x <- Value_raw[,i]
    valid_index <- which(x>0)
    valid_x <- x[valid_index]
    valid_x_sum <- sum(valid_x)
    valid_x_FOT5 <- valid_x/valid_x_sum*1e5
    Value_FOT5[valid_index,i] <- valid_x_FOT5
  }
  data_frame_normaliation <- data.frame(ID, Value_FOT5)
  data_frame_normaliation_colnames <- c(data_frame_colnames[1], experiment_code)
  colnames(data_frame_normaliation) <- data_frame_normaliation_colnames
  return(data_frame_normaliation)
}
