#' Motif enrichment using rmotifx.
#'
#' @param foreground A vector for aligned sequences as the foreground input.
#' @param background A vector for aligned sequences as the background input.
#' @param center_vector A vector for aligned centers.
#' @param motifx_pvalue A numeric value for selecting motifs that meets the minimum cutoff.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Omar Wagih (2014). rmotifx: An iterative statistical approach to the discovery of biological sequence motifs. R package version 1.0.
#'
#' @return A list for results of motif enrichment.
#' @export
#'
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#'     ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_motifs_list.RData"
#'     load_data <- load_data_with_ftp(ftp_url, 'RData')
#'     writeBin(load_data, "get_motifs_list.RData")
#'     load("get_motifs_list.RData")
#'
#'     motifs_list <- get_motifs_list(foreground[1:100], background[1:100], center_vector, motifx_pvalue)
#'     head(motifs_list)
#' }


get_motifs_list <- function(foreground, background, center_vector, motifx_pvalue){
  motifs_list <- list()
  motifs_list_names <- NULL
  motifs_list_index <- 0
  center_vector_len <- length(center_vector)
  for(i in seq_len(center_vector_len)){
    center <- center_vector[i]
    motifs <- get_motif_analysis_summary(foreground, background, center = center, min_sequence_count = 1, min_pvalue = motifx_pvalue)
    if(!is.null(motifs)){
      motifs_list_index <- motifs_list_index + 1
      motifs_list[[motifs_list_index]] <- motifs
      motifs_list_names <- c(motifs_list_names, center)
    }
  }
  if(motifs_list_index > 0){
    names(motifs_list) <- motifs_list_names
    return(motifs_list)
  }else{
    return(NULL)
  }
}
