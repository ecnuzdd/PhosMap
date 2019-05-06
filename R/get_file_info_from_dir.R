#' Get data lists from and corresponding file ids.
#'
#' Read batch files (.txt or .csv) from specific directory.
#'
#' @param specific_dir A folder containing files as input.
#' @param experiment_ID A vector containing experiment codes as input
#'
#' @return A list containing data from files and corresponding file ids
#' @export
#'
#' @examples
#' \dontrun{
#' result_list = get_file_info_from_dir(
#'   specific_dir,
#'   experiment_ID
#' )
#' }
get_file_info_from_dir <- function(specific_dir, experiment_ID){
  requireNamespace('utils')
  # read all files from specific director and save them into a list
  all_files <- list.files(specific_dir)
  all_files_count <- length(all_files)
  if(all_files_count>0){
    file_suffix <- get_file_suffix(all_files[1])
    if(file_suffix=='txt'){
      read_file_function <- utils::read.table
      sep <- '\t'
    }else{
      read_file_function <- utils::read.csv
      sep <- ','
    }
    sep_symbol <- paste('.', file_suffix, sep = '')
    all_files_ID <- apply(data.frame(all_files), 1, function(x, sep){
      x <- strsplit(x, split = sep)[[1]][1]
      x
    }, sep=sep_symbol)

    all_files_ID_code <- apply(data.frame(all_files_ID), 1, function(x, sep){
      x <- strsplit(x, split = sep)[[1]][1]
      x
    }, sep='_')
    all_files_paths <- normalizePath(file.path(specific_dir, all_files))

    index_of_match <- match(experiment_ID, all_files_ID_code)
    matched_all_files_paths <- all_files_paths[index_of_match]
    matched_all_files_ID <- all_files_ID[index_of_match]

    file_data_list <- list()
    matched_all_files_count <- length(matched_all_files_paths)
    cat('\n Total file: ', matched_all_files_count)
    for(i in seq_len(matched_all_files_count)){
      # Read bach data and save to file_data_list.
      cat('\n completed: ', i, '/', matched_all_files_count)
      file_data <- as.matrix(read_file_function(matched_all_files_paths[i], header = TRUE, sep = sep))
      file_data_list[[i]] <- file_data
    }
    attr(file_data_list,'names') <- matched_all_files_ID
    result_list <- list(file_data_list=file_data_list, file_ID=matched_all_files_ID)
    return(result_list)

  }else{
    stop('The directory of ', specific_dir, ' has no files.')
  }
}


