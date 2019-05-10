#' Get background data frame (fasta library from Refseq).
#'
#' @param species A string for that the alignment is based on which species, the options are human, mouse and rattus.
#' @param fasta_type, A string for fasta source, the options are refseq and uniprot, the default is refseq
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame of background
#' @export
#'
#' @examples
#' \dontrun{
#' background_df <- get_global_background_df(species)
#' }


get_global_background_df <- function(species = 'human', fasta_type = 'refseq'){
  requireNamespace('utils')

  # *** background ***
  cat('Reading background.\n')
  if(fasta_type == 'refseq' | fasta_type == 'uniprot'){
    # Read motif background
    BACKGROUND_FILE_DIR <- normalizePath(
      system.file(
        'extdata',
        'motif_library', fasta_type, species,
        package = "PhosMap"
      ),
      mustWork = F
    )

  }else{
    cat('\n The fasta type is incorrect, the expected type is refseq or uniprot.')
    stop('')
  }

  background_file_name <- paste('STY_background_of_',fasta_type, '_', species, '_for_motif_enrichment.txt', sep = '')
  BACKGROUND_FILE_PATH <- normalizePath(
    file.path(BACKGROUND_FILE_DIR, background_file_name), mustWork = FALSE
  )
  if(!file.exists(BACKGROUND_FILE_PATH)){
    PHOSPHATE_LIB_MOTIF_BG_FILE_http_link <- 'https://media.githubusercontent.com/media/ecnuzdd/PhosMap_datasets/master/motif_library/fasta_type/species/STY_background_of_fasta_type_species_for_motif_enrichment.txt'
    PHOSPHATE_LIB_MOTIF_BG_FILE_http_link <- stringr::str_replace_all(PHOSPHATE_LIB_MOTIF_BG_FILE_http_link, 'fasta_type', fasta_type)
    PHOSPHATE_LIB_MOTIF_BG_FILE_http_link <- stringr::str_replace_all(PHOSPHATE_LIB_MOTIF_BG_FILE_http_link, 'species', species)

    PHOSPHATE_LIB_MOTIF_BG_FILE_data_type <- 'txt'

    PHOSPHATE_LIB_MOTIF_BG_DATA <- tryCatch(
      {
        load_data_with_http(PHOSPHATE_LIB_MOTIF_BG_FILE_http_link, PHOSPHATE_LIB_MOTIF_BG_FILE_data_type)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        return('error')
      },
      warning=function(cond) {
        message("Here's the original warning message:")
        message(cond)
        return('warning')
      }
      # ,
      # finally={
      #   message("Some other message at the end")
      # }
    )

    if(is.list((PHOSPHATE_LIB_MOTIF_BG_DATA))){
      cat('\n')
      message('The current network is unstable when downloading background data from ', PHOSPHATE_LIB_MOTIF_BG_FILE_http_link)
      cat('\n')
      message('You can try it again.')
      cat('\n')
      message('The another option is that you can directly download the background data from ', PHOSPHATE_LIB_MOTIF_BG_FILE_http_link, ' and then put it into ', BACKGROUND_FILE_PATH)
      stop()
    }
    message('Save ', fasta_type, ' fasta file of ', species, ' to ', BACKGROUND_FILE_PATH)
    write.table(PHOSPHATE_LIB_MOTIF_BG_DATA, BACKGROUND_FILE_PATH, row.names = FALSE, col.names = TRUE, sep = '\t')
    message('Save successfully.')
  }else{
    cat('Read background file of ', species, '.\n', sep = '')
    PHOSPHATE_LIB_MOTIF_BG_DATA <- utils::read.table(BACKGROUND_FILE_PATH, sep = '\t', header = TRUE)
    cat('Read OK! ^_^', '\n')
  }
  background_df <- PHOSPHATE_LIB_MOTIF_BG_DATA
  return(background_df)

}
