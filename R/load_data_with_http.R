#' Load datasets from URL (https://github.com/ecnuzdd/PhosMap_datasets)
#'
#' Some datasets with larger size need to be loaded for mapping ID or protein sequence when using PhosMap.
#' These datasets could be ragarded as library and uploaded to https://github.com/ecnuzdd/PhosMap_datasets in advance.
#' When first perfoming functions depending on these datasets, PhosMap will get them from specific URL and save them into local disk.
#'
#'
#' @param http_link A string for URL of datasets.
#' @param data_type A string for type of datasets (txt of csv).
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' load_data <- load_data_with_http(
#'   http_link, data_type
#' )
#' }

load_data_with_http <- function(http_link, data_type){
  message('First loading data from Github sever, it may take a few minutes.')
  message('Downloading data from ', http_link, '.', sep = '')
  # message('Suggesting you save the data for next use.')
  if(data_type == 'csv'){
    load_data = read.csv(http_link)
    # save
  }else{
    load_data = read.table(http_link, header = T, sep = '\t')
    # save
  }
  message('Completing the data load.')
  return(load_data)
}
