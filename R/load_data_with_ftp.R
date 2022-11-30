#' Load datasets from URL (ftp://111.198.139.72:4000/pub/PhosMap_datasets)
#'
#' Some datasets with larger size need to be loaded for mapping ID or protein sequence when using PhosMap.
#' These datasets could be ragarded as library and uploaded to ftp://111.198.139.72:4000/pub/PhosMap_datasets in advance.
#' When first perfoming functions depending on these datasets, PhosMap will get them from specific URL and save them into local disk.
#'
#'
#' @param ftp_link A string for URL of datasets.
#' @param data_type A string for type of datasets (txt, csv, RData).
#' @import RCurl utils
#' @author Dongdong Zhan and Mengsha Tong
#' 
#' @return A dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' ftp_url <- "https://github.com/ecnuzdd/PhosMap_datasets/function_demo_data/profiling_exp_design_info.txt"
#' load_data <- load_data_with_ftp(ftp_url, 'txt')
#' head(load_data)
#' }

load_data_with_ftp <- function(ftp_link, data_type){
  # ftp_url <- "https://github.com/ecnuzdd/PhosMap_datasets/BRAFi.RData"
  # load_data <- load_data_with_ftp(ftp_url, 'RData')
  # writeBin(load_data, "BRAFi.RData")
  # load("BRAFi.RData")

  # ftp_url <- "https://github.com/ecnuzdd/PhosMap_datasets/motif_library/refseq/mouse/STY_background_of_refseq_mouse_for_motif_enrichment.txt"
  # bg <- load_data_with_ftp(ftp_url, 'txt')

  # ftp_url <- "https://github.com/ecnuzdd/PhosMap_datasets/kinase_substrate_regulation_relationship_table/human/human_ksrr.csv"
  # ks <- load_data_with_ftp(ftp_url, 'csv')
  requireNamespace('RCurl')
  ftp_url <- ftp_link
  message('First loading data from FTP sever, it may take a few minutes.')
  message('Downloading data from ', ftp_link, '.', sep = '')
  # message('Suggesting you save the data for next use.')
  userpwd <- "user1:qazwsx123!"
  if(data_type == 'csv'){
    csv_data <- getURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE, crlf = TRUE)
    load_data = read.csv(text = csv_data)
    message('Completing the csv data load.')
    return(load_data)
  }else if(data_type == 'txt'){
    text_data <- getURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE, crlf = TRUE)
    load_data <- read.table(text = text_data, header = TRUE, sep = '\t')
    message('Completing the text data load.')

  }else{
    load_data <- RCurl::getBinaryURL(ftp_url, userpwd = userpwd, ftp.use.epsv = FALSE, crlf = TRUE)
    message('Completing the RData load.')
  }
  return(load_data)
}
