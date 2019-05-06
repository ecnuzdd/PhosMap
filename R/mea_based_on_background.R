#' Motif enrichment based on global background (fasta library from Refseq).
#'
#' @param foreground A vector for aligned sequence of foreground.
#' @param AA_in_protein A vector for the location of S/T/Y in sequence of protein.
#' @param background A vector for aligned sequence of background.
#' @param motifx_pvalue A numeric value for selecting motifs that meets the minimum cutoff.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A list containing motifs and the corresponding sequences
#' @export
#'
#' @examples
#' \dontrun{
#' motifs_list = mea_based_on_background(
#'   foreground,
#'   AA_in_protein,
#'   background,
#'   motifx_pvalue
#' )
#' }



# foreground = as.vector(foreground_df$aligned_seq)
# AA_in_protein = as.vector(foreground_df$AA_in_protein)
mea_based_on_background <- function(foreground, AA_in_protein, background, motifx_pvalue){
  # foreground = as.vector(foreground)
  # background = as.vector(background$Aligned_Seq)
  center_vector_candidate <- c('S', 'T', 'Y')
  center_vector_candidate_len <- length(center_vector_candidate)
  center_vector <- NULL
  for(i in seq_len(center_vector_candidate_len)){
    center <- center_vector_candidate[i]
    if(length(grep(center, AA_in_protein)) > 0){
      center_vector <- c(center_vector, center)
    }
  }
  cat('Start executing motifx and find motif pattern. \n')
  cat('Foreground sequences: ', length(foreground), '.\n', sep = '')
  cat('Background sequences: ', length(background), '.\n', sep = '')
  cat('Phosphorylation: [', center_vector, '] exists in foreground.\n', sep = '')
  cat('Motifx pvalue cutoff: ', motifx_pvalue, '.\n', sep = '')
  motifs_list <- get_motifs_list(foreground, background, center_vector, motifx_pvalue)
  cat('Motifx analysis OK! ^_^', '\n')
  print(motifs_list)
  cat('\n')
  return(motifs_list)
}













































