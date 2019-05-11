#' computing kinase-substrate enrichment score
#'
#' @param substate_vector a vector for substrates with values indentified in current experiments.
#' @param regulons_of_kinase a vector for substrates of a specific kinase, which with substrates identified in current experiments.
#' @param substrates_of_kinase_in_exp_count a numeric for numbers in regulons_of_kinase vector.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Hernandez-Armenta C et al. Benchmarking substrate-based kinase activity inference using phosphoproteomic data[J]. Bioinformatics, 2017, 33(12): 1845-1851.
#'
#' @return A numeric or NA for enrichment_score.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' enrichment_score <- compute_kses(
#'   substate_vector,
#'   regulons_of_kinase,
#'   substrates_of_kinase_in_exp_count
#' )
#' }
#'

compute_kses <- function(
  substate_vector,
  regulons_of_kinase,
  substrates_of_kinase_in_exp_count
){
  if(substrates_of_kinase_in_exp_count > 0){
    substrates_value_in_exp <- substate_vector
    substrates_name_in_exp <- names(substate_vector)
    match_flag <- is.element(substrates_name_in_exp, regulons_of_kinase)
    substrates_name_in_exp_count <- length(substrates_name_in_exp)

    match_profile <- rep(0, substrates_name_in_exp_count)
    names(match_profile) <- substrates_name_in_exp
    match_profile[match_flag] <- substrates_value_in_exp[match_flag]
    match_profile_cumsum <- cumsum(abs(match_profile))
    match_profile_cumsum_max <- max(match_profile_cumsum)

    unmatch_profile <- rep(1, substrates_name_in_exp_count)
    names(unmatch_profile) <- substrates_name_in_exp
    unmatch_profile[match_flag] <- 0
    unmatch_profile_cumsum <- cumsum(abs(unmatch_profile))

    non_substrates_of_kinase_in_exp_count <- substrates_name_in_exp_count - substrates_of_kinase_in_exp_count

    match_probs <- match_profile_cumsum/match_profile_cumsum_max
    unmatch_probs <- unmatch_profile_cumsum/non_substrates_of_kinase_in_exp_count
    probs_diff <- abs(match_probs - unmatch_probs)

    probs_diff_max_index <- which.max(probs_diff)[1]

    enrichment_score <- match_probs[probs_diff_max_index] - unmatch_probs[probs_diff_max_index]

  }else{
    enrichment_score <- NA
  }
  return(enrichment_score)
}
