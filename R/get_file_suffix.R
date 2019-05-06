#' Get sufffix of input file.
#'
#' @param file_name A string for file name.
#'
#' @return A string for file format.
#' @export
#'
get_file_suffix <- function(file_name){
  # Paser file name and get file suffix.
  # txt or csv
  if(grepl('.txt', file_name)){
    return('txt')
  }else{
    return('csv')
  }
}
