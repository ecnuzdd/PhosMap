#' Convert data frame of motif to the sequence pattern
#'
#' @param motif_data_frame A data frame with two columns including amino acid and index on sequence with fixed length.
#' @param center A character for center of k-mer.
#' @param width A numeric for specific k-mer.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#' @return A string for motif pattern
#' @export
#'
#' @examples
#' \dontrun{
#' motif_pattern <- motif_data_frame_to_sequence(
#'   motif_data_frame,
#'   center,
#'   width
#' )
#' }

motif_data_frame_to_sequence <- function(motif_data_frame, center, width){
  motif_pattern <- rep('.', width)
  motif_pattern[ceiling(width/2)] <- center
  motif_data_frame_row <- nrow(motif_data_frame)
  for(i in seq_len(motif_data_frame_row)){
    motif_pattern[motif_data_frame[i,2]] <- motif_data_frame[i,1]
  }
  motif_pattern <- paste(motif_pattern, collapse = '')
  return(motif_pattern)
}
