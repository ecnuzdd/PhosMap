#' Get indexes of modifications in protein sequences.
#' @param id_data_only_peptide2gi a data frame for peptides with protein gi.
#' @param fasta_data a fasta data for a specific species.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A vector for indexes of modifications in protein sequences.
#' @export
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_modification_index.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_modification_index.RData")
#' load("get_modification_index.RData")
#'
#' modification_index_in_protein_seq_list <- get_modification_index(
#'   id_data_only_peptide2gi[1:100, ],
#'   fasta_data
#' )
#' head(modification_index_in_protein_seq_list)
#' }
#'

get_modification_index <- function(id_data_only_peptide2gi, fasta_data){
  # 1
  # Get modification index in protein sequence.
  cat('\n', 'Get modification index in protein sequence.')
  id_data_only_peptide2gi_row <- nrow(id_data_only_peptide2gi)
  modification_index_in_protein_seq_list <- list()
  for(i in seq_len(id_data_only_peptide2gi_row)){
    peptide_seq <- as.vector(id_data_only_peptide2gi$Sequence[i])
    peptide_id <- as.vector(id_data_only_peptide2gi$ID[i])
    modification_index_in_peptide_seq <- unlist(gregexpr("[a-z]", peptide_seq))
    protein_seq <- as.vector(fasta_data$Sequence[which(fasta_data$ID==peptide_id)])
    first_index_of_peptide2protein <- unlist(gregexpr(toupper(peptide_seq), protein_seq))
    modification_index_in_protein_seq <- NULL
    for(elemt in first_index_of_peptide2protein){
      tmp_modification_index_in_protein_seq <- elemt + modification_index_in_peptide_seq -1
      modification_index_in_protein_seq <- c(modification_index_in_protein_seq,
                                             tmp_modification_index_in_protein_seq)
    }
    modification_index_in_protein_seq_list[[i]] <- modification_index_in_protein_seq
    if(i%%500==0 | i==id_data_only_peptide2gi_row ){
      cat('\n completed: ', i, '/', id_data_only_peptide2gi_row)
    }
  }
  return(modification_index_in_protein_seq_list)
}
