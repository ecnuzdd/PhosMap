#' Load datasets from URL (ftp://111.198.139.72:4000/pub/PhosMap_datasets)
#'
#' Some datasets with larger size need to be loaded for mapping ID or protein sequence when using PhosMap.
#' These datasets could be ragarded as library and uploaded to ftp://111.198.139.72:4000/pub/PhosMap_datasets in advance.
#' When first perfoming functions depending on these datasets, PhosMap will get them from specific URL and save them into local disk.
#'
#'
#' @param ftp_link A string for URL of datasets.
#' @param data_type A string for type of datasets (txt, csv, RData).
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' http_link <- url('https://raw.githubusercontent.com/ecnuzdd/PhosMap/master/inst/extdata/kinase_substrate_regulation_relationship_table/human/human_ksrr.csv')
#' data_type = 'csv'
#' load_data <- load_data_with_http(
#'   http_link, data_type
#' )
#'

load_data_with_ftp <- function(ftp_link, data_type){
  message('First loading data from Github sever, it may take a few minutes.')
  message('Downloading data from ', http_link, '.', sep = '')
  # message('Suggesting you save the data for next use.')

  require('RCurl')
  userpwd <- "user1:qazwsx123!"
  if(data_type == 'csv'){
    csv_data <- getURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE, crlf = TRUE)
    load_data = read.csv(csv_data)
    message('Completing the csv data load.')
    return(load_data)
  }else if(txt){
    text_data <- getURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE, crlf = TRUE)
    load_data = read.table(text_data, header = TRUE, sep = '\t')
    message('Completing the text data load.')
    return(load_data)
  }else{
    out <- getBinaryURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE,crlf = TRUE)
    writeBin(out, "temp.RData")
    load("temp.RData")
    message('Completing the RData load.')
  }
}
