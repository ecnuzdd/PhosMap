#' computing kinase-substrate enrichment significance (pvalue)
#'
#' @param substate_vector a vector for substrates with value indentified in current experiment.
#' @param regulons_of_kinase a vector for substrates of a specific kinase, which identified in current experiment.
#' @param trial a numeric for the number of random samples, the default is 1000.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Hernandez-Armenta C et al. Benchmarking substrate-based kinase activity inference using phosphoproteomic data[J]. Bioinformatics, 2017, 33(12): 1845-1851.
#'
#' @return A list for expected enrichment score and its significance
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pvalue <- get_kses(
#'   substate_vector,
#'   regulons_of_kinase,
#'   1000
#' )
#' }
#'
#'
get_kses <- function(
  substate_vector,
  regulons_of_kinase,
  trial = 1000
){
  substrates_of_kinase_in_exp <- intersect(names(substate_vector), regulons_of_kinase)
  substrates_of_kinase_in_exp_count <- length(substrates_of_kinase_in_exp)
  expected_enrichment_score <- compute_kses(
    substate_vector,
    regulons_of_kinase,
    substrates_of_kinase_in_exp_count
  )

  if(is.na(as.vector(expected_enrichment_score))){
    return(list())
  }

  stochastic_enrichment_scores <- NULL
  substate_vector_count <- length(substate_vector)
  # stochastic process
  for(i in seq_len(trial)){
    regulons_of_kinase_i <- names(substate_vector[sample(substate_vector_count, substrates_of_kinase_in_exp_count)])
    # print(regulons_of_kinase_i)
    stochastic_enrichment_score_i <- compute_kses(
      substate_vector,
      regulons_of_kinase_i,
      substrates_of_kinase_in_exp_count
    )
    stochastic_enrichment_scores <- c(stochastic_enrichment_scores, stochastic_enrichment_score_i)
  }
  if(as.vector(expected_enrichment_score) >= 0){
    pvalue <- length(which(stochastic_enrichment_scores>=expected_enrichment_score))/trial
  }else{
    pvalue <- length(which(stochastic_enrichment_scores<expected_enrichment_score))/trial
  }
  return(list(pvalue = pvalue, expected_enrichment_score = expected_enrichment_score))
}
