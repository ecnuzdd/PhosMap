#' Check input for motif enrichment analysis (mea)
#'
#' @param foreground A vector for AA sequences with fixed length as foreground input.
#' @param background A vector for AA sequences with fixed length as background input.
#' @param center A character for center of k-mer.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#' @return A list passing check steps
#' @export
#'
#' @examples
#' \dontrun{
#' check_result_list <- check_mea_input(
#'   foreground,
#'   background,
#'   center
#' )
#' }

check_mea_input <- function(
  foreground,
  background,
  center
){
  AA_LIST <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

  if(length(foreground) == 0){
    stop('No foreground sequences.')
  }
  if(length(background) == 0){
    stop('No background sequences.')
  }

  if(!(center %in% AA_LIST)){
    stop('The center is not an amino acid.')
  }

  width = nchar(foreground[1])
  if(width %% 2 == 0){
    stop('Sequence length must be an odd number.')
  }
  if(width < 3 | width > 25){
    stop('Sequence length is limited to between 3 and 25.')
  }
  if(any(nchar(foreground) != width)){
    stop('All foreground must be the same length.')
  }
  if(width != nchar(background[1])){
    stop('Background length must be same as foreground.')
  }
  if(any(nchar(background) != width)){
    stop('All background must be the same length.')
  }
  center_index <- ceiling(width/2)
  foreground <- foreground[substr(foreground, center_index, center_index) == center]
  background <- background[substr(background, center_index, center_index) == center]
  if(length(foreground)>0){
    foreground <- foreground[!grepl('[BJOUXZ]', foreground)]
  }
  if(length(background)>0){
    background <- background[!grepl('[BJOUXZ]', background)]
    background <- unique(background)
  }
  check_result_list <- list(
    foreground = foreground,
    background = background,
    width = width
  )
}
