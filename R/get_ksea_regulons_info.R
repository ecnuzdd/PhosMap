#' Get informational data frame by combining results from all experiments
#'
#' @param ksea_regulons A kinase vector from all experiments.
#' @param ksea_trans_list A list that consits of regulation direction of kinase from each experiment by ksea.
#' @param ksea_x_list A list that consits of sepecific information from each experiment by ksea, like regulation direction, p value and activity etc..
#' @param ptypes_data_ratio_colnames A vector that consists of column names from experiments.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame containing sepecific information of all experiments from ksea results, like regulation direction, pvalue and activity etc..
#'
#' @export
#'
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_ksea_regulons_info.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_ksea_regulons_info.RData")
#' load("get_ksea_regulons_info.RData")
#'
#' ksea_regulons_activity_df <- get_ksea_regulons_info(
#'   ksea_regulons,
#'   ksea_trans_list,
#'   ksea_activity_list,
#'   ptypes_data_ratio_colnames
#' )
#' ksea_regulons_activity_df
#'

get_ksea_regulons_info <- function(
  ksea_regulons,
  ksea_trans_list,
  ksea_x_list,
  ptypes_data_ratio_colnames
){
  ksea_regulons_count <- length(ksea_regulons)
  ksea_regulons_x_df <-  c()
  ptypes_data_exp_count <- length(ptypes_data_ratio_colnames)
  for(i in seq_len(ptypes_data_exp_count)){
    ksea_trans_i <- ksea_trans_list[[i]]
    ksea_trans_i_names <- names(ksea_trans_i)
    ksea_x_list_i <- ksea_x_list[[i]]
    ksea_regulons_x_i <- apply(data.frame(ksea_regulons), 1, function(x, ksea_x_list_i, ksea_trans_i_names){
      index_of_match <- which(ksea_trans_i_names==x)
      if(length(index_of_match)>0){
        return(ksea_x_list_i[index_of_match])
      }else{
        return(0)
      }
    }, ksea_x_list_i = ksea_x_list_i, ksea_trans_i_names = ksea_trans_i_names)
    ksea_regulons_x_df <- rbind(ksea_regulons_x_df, ksea_regulons_x_i)
  }
  colnames(ksea_regulons_x_df) <- ksea_regulons
  rownames(ksea_regulons_x_df) <- ptypes_data_ratio_colnames
  Kinase <- ksea_regulons
  ksea_regulons_x_df <- t(ksea_regulons_x_df)
  ksea_regulons_x_df <- data.frame(Kinase, ksea_regulons_x_df)
  return(ksea_regulons_x_df)
}
