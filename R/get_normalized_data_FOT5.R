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
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url1 <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_normalized_data_FOT5.RData"
#' ftp_url2 <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/profiling_exp_design_info.txt"

#' load_data1 <- load_data_with_ftp(ftp_url1, 'Rdata')
#' writeBin(load_data1, "get_normalized_data_FOT5.RData")
#' load("get_normalized_data_FOT5.RData")
#'
#' load_data2 <- load_data_with_ftp(ftp_url2, 'downloadtxt')
#' writeBin(load_data2, "profiling_exp_design_info.txt")
#' profiling_exp_design_info_file_path <- "./profiling_exp_design_info.txt"
#'
#' profiling_data_normalized <- get_normalized_data_FOT5(profiling_data,
#'   profiling_exp_design_info_file_path
#' )
#' head(profiling_data_normalized)
#'
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
