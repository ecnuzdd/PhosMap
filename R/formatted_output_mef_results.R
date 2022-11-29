#' Output formatted sequences in foreground that are mapped to specific motifs.
#'
#' @param foreground_sequences_mapped_to_motifs A list that consists of motifs and their corresponding aligned sequences from foreground.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#' @return A data frame that motifs and their corresponding aligned sequences from foreground.
#' @export
#'
#' @examples
#' \dontrun{
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/formatted_output_mef_results.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "formatted_output_mef_results.RData")
#' load("formatted_output_mef_results.RData")
#'
#' formatted_output_df <- formatted_output_mef_results(
#'   foreground_sequences_mapped_to_motifs
#' )
#' head(formatted_output_df)
#' }

formatted_output_mef_results <- function(
  foreground_sequences_mapped_to_motifs
){
  cat('Output formatted sequences in foreground that are mapped to specific motif. \n')
  formatted_output_vector <- NULL
  foreground_sequences_mapped_to_motifs_len <- length(foreground_sequences_mapped_to_motifs)
  for(i in seq_len(foreground_sequences_mapped_to_motifs_len)){
    motif <- names(foreground_sequences_mapped_to_motifs[i])
    mapped_seq <- foreground_sequences_mapped_to_motifs[[i]]
    formatted_output_vector <- c(formatted_output_vector, motif, mapped_seq, '')
  }
  formatted_output_df <- data.frame(formatted_output_vector)
  cat('\n')
  return(formatted_output_df)
}
