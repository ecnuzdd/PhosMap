#' Get summary results of motif analysis for specific input
#'
#' @param foreground A vector for AA sequences with fixed length as foreground input.
#' @param background A vector for AA sequences with fixed length as background input.
#' @param center A character for center of k-mer.
#' @param min_sequence_count A numeric for the minimum sequence number assigned to a motif.
#' @param min_pvalue A numeric for the minimum pvalue for found motif.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @references Omar Wagih (2014). rmotifx: An iterative statistical approach to the discovery of biological sequence motifs. R package version 1.0.
#'
#'
#' @return A list for summary result of motif analysis
#' @export
#'
#' @examples
#' \dontrun{
#' summry_list <- get_motif_analysis_summary(
#'   foreground,
#'   background,
#'   center,
#'   min_sequence_count,
#'   min_pvalue
#' )
#' }

get_motif_analysis_summary <- function(
  foreground,
  background,
  center='S',
  min_sequence_count = 1,
  min_pvalue = 0.01
){
  check_result_list <- check_mea_input(foreground, background, center)
  loop_foreground <- check_result_list$foreground
  loop_background <- check_result_list$background
  motif_result_list <- list()
  motif_result_list_index <- 0
  while(length(loop_foreground) >= min_sequence_count){
    motif_result_loop_i <- seach_motif_pattern(
      loop_foreground,
      loop_background,
      min_sequence_count = min_sequence_count,
      min_pvalue = min_pvalue,
      center = center,
      width = check_result_list$width
    )
    if(is.null(motif_result_loop_i)){
      break
    }
    motif_result_list_index <- motif_result_list_index + 1
    motif_result_list[[motif_result_list_index]] <- motif_result_loop_i
    loop_foreground <- loop_foreground[!grepl(motif_result_loop_i$motif_pattern, loop_foreground)]
    loop_background <- loop_background[!grepl(motif_result_loop_i$motif_pattern, loop_background)]
  }

  summry_list <- data.frame(
    motif = vapply(motif_result_list, function(x){x$motif_pattern},c('character')),
    score = vapply(motif_result_list, function(x){x$motif_pattern_score}, c(1)),
    foreground_matches = vapply(motif_result_list, function(x){x$foreground_matches}, 1),
    foreground_size = vapply(motif_result_list, function(x){x$foreground_size}, 1),
    background_matches = vapply(motif_result_list, function(x){x$background_matches}, 1),
    background_size = vapply(motif_result_list, function(x){x$background_size}, 1)
  )

  foreground_fold_increase <- summry_list$foreground_matches/summry_list$foreground_size
  background_fold_increase <- summry_list$background_matches/summry_list$background_size
  summry_list$fold_increase <- foreground_fold_increase/background_fold_increase

  if(nrow(summry_list) == 0){
    return(NULL)
  }
  return(summry_list)
}
