#' Keep psites whose row maximum is top N (percentage).
#'
#' Compute row maximum each psites, sort row maximum in decreasing order and keep top N (percentage).
#'
#' @param phospho_data A data frame of phospho-data.
#' @param percent_of_kept_sites A numeric value representing a cutoff used for filter psites. The default is 3/4.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame meeting specific cutoff.
#' @export
#'
#' @examples
#' \dontrun{
#' phospho_data_meet_percent = keep_psites_with_max_in_topX(
#'   phospho_data,
#'  percent_of_kept_sites = 3/4
#' )
#' }
keep_psites_with_max_in_topX <- function(phospho_data, percent_of_kept_sites = 3/4){
  percent_of_kept_sites_str <- paste('top', percent_of_kept_sites*100, '%', sep = '')
  cat('\n The 8th step: filter psites with row maximum in', percent_of_kept_sites_str, '.')
  ID <- as.vector(phospho_data[,1])
  Value <- phospho_data[,-1]
  Value_rowmax <- apply(Value, 1, function(x){
    x <- as.vector(unlist(x))
    max(x)
  })
  index_of_Value_rowmax_desc <- order(Value_rowmax, decreasing = TRUE)
  count_of_kept_sites <- round(nrow(Value)*percent_of_kept_sites)
  index_of_Value_rowmax_desc_kept <- index_of_Value_rowmax_desc[seq_len(count_of_kept_sites)]
  phospho_data_meet_percent <- phospho_data[index_of_Value_rowmax_desc_kept,]
  cat('\n The 8th step: filter over with ', percent_of_kept_sites_str, ' cutoff.')
  return(phospho_data_meet_percent)
}
