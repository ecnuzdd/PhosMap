#' Output formatted sequences in foreground that are mapped to specific motif.
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
#' formatted_output_df <- formatted_output_mef_results(
#'   foreground_sequences_mapped_to_motifs
#' )
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
