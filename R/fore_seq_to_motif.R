#' Convert the list that consists of motifs and the corresponding sequences to data frame.
#'
#' @param foreground_sequences_mapped_to_motifs A list that consists of motifs and the corresponding sequences.
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame that consist of aligned sequences and the corresponding motifs.
#' @export
#'
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/fore_seq_to_motif.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "fore_seq_to_motif.RData")
#' load("fore_seq_to_motif.RData")
#'
#' df <- fore_seq_to_motif(
#'   foreground_sequences_mapped_to_motifs
#' )
#' head(df)
#'

fore_seq_to_motif <- function(
  foreground_sequences_mapped_to_motifs
){
  foreground_sequences_mapped_to_motifs_len <- length(foreground_sequences_mapped_to_motifs)
  motif_rep_v <- NULL
  map_seq_v <- NULL
  for(i in seq_len(foreground_sequences_mapped_to_motifs_len)){
    motif <- names(foreground_sequences_mapped_to_motifs[i])
    map_seq <- foreground_sequences_mapped_to_motifs[[i]]
    motif_rep <- rep(motif, length(map_seq))
    motif_rep_v <- c(motif_rep_v, motif_rep)
    map_seq_v <- c(map_seq_v, map_seq)
  }
  df <- data.frame(map_seq_v, motif_rep_v)
  colnames(df) <- c('Aligned_seq', 'Motif')
  return(df)
}
