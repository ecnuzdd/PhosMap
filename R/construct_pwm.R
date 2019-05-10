#' Construct position weight matrix
#'
#' @param sequences A vector for aligned sequences with fixed length.
#' @param width A numeric for specific k-mer.
#' @param frequency_flag A boolean for showing real frequency or frequency probability, the default is TRUE for showing real frequency.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Wagih O, Sugiyama N, Ishihama Y, et al. Uncovering phosphorylation-based specificities through functional interaction networks[J]. Molecular & Cellular Proteomics, 2016, 15(1): 236-245.
#'
#' @return A position weight matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' pwm <- construct_pwm(
#'   sequences,
#'   width,
#'   frequency_flag = TRUE
#' )
#' }

construct_pwm <- function(
  sequences,
  width,
  frequency_flag = TRUE
){
  AA_LIST <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  pwm <- matrix(0, length(AA_LIST), width, dimnames = list(AA_LIST, seq(1, width)))
  sequences_count <- length(sequences)
  for(i in seq_len(sequences_count)){
    sequences_i <- sequences[i]
    for(j in seq_len(width)){
      char_j <- substr(sequences_i, j, j)
      if(char_j != "_"){
        pwm[char_j, j] <- pwm[char_j, j] + 1
      }
    }
  }
  if(!frequency_flag){
    pwm_colsums <- colSums(pwm)
    for(i in seq_len(width)){
      pwm[,i] <- pwm[,i]/pwm_colsums[i]
    }
  }
  return(pwm)
}
