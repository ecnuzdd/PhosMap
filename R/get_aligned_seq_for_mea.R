#' Taking S/T/Y as the center, align sequence to fasta library by the specific length.
#'
#' @param ID A vector for gi number of proteins.
#' @param Sequence A vector for sequence of peptides.
#' @param AA_in_protein A vector for the locations of S/T/Y in sequence of proteins.
#' @param fixed_length Length of aligned sequence,the default is 15.
#' @param species A string for the library of species, the options are human, mouse and rattus, the default is human.
#' @param fasta_type, A string for fasta source, the options are refseq and uniprot, the default is refseq
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.3.0.\
#' https://CRAN.R-project.org/package=stringr.
#' @import utils
#' @return A data frame containing ID, Sequence, AA_in_protein, aligned_seq.
#' @export
#'
#' @examples
#' ## The process needs to load data from PhosMap datasets stored into FTP server and perform large computation.
#' ## It may take a few minutes.
#' if(FALSE){
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_aligned_seq_for_mea.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_aligned_seq_for_mea.RData")
#' load("get_aligned_seq_for_mea.RData")
#'
#' foreground_df <- get_aligned_seq_for_mea(
#'   ID[1:100], Sequence[1:100], AA_in_protein[1:100],
#'   fixed_length, species = 'human',
#'   fasta_type = 'refseq'
#' )
#' head(foreground_df)
#'
#' }



get_aligned_seq_for_mea <- function(ID, Sequence, AA_in_protein, fixed_length, species = 'human', fasta_type = 'refseq'){
  requireNamespace('stringr')
  requireNamespace('utils')
  # require(PhosMap)
  cat('Aligned sequence based on fasta library for motif enrichment anlysis.\n')

  ###########################################################################################################
  # Read fasta file: fasta_data
  cat('Read fasta file of ', species, '.\n', sep = '')
  # fasta_type <- 'refseq' # 'uniprot'
  if(fasta_type == 'refseq' | fasta_type == 'uniprot'){
    PHOSPHATE_LIB_FASTA_DIR <- normalizePath(
      system.file(
        'extdata', 'fasta_library',
        fasta_type, species,
        package = "PhosMap"
      ),
      mustWork = FALSE
    )

    PHOSPHATE_LIB_FASTA_FILE_PATH <- normalizePath(
      file.path(PHOSPHATE_LIB_FASTA_DIR, paste(species, fasta_type, 'fasta.txt', sep = '_')),
      mustWork = FALSE
    )
  }else{
    cat('\n The fasta type is incorrect, the expected type is refseq or uniprot.')
    stop('')
  }

  # 'ftp://111.198.139.72:4000/pub/PhosMap_datasets/fasta_library/refseq/human/human_ref_fasta.txt
  if(!file.exists(PHOSPHATE_LIB_FASTA_FILE_PATH)){
    PHOSPHATE_LIB_FASTA_FILE_ftp_link <- 'ftp://111.198.139.72:4000/pub/PhosMap_datasets/fasta_library/fasta_type/species/species_fasta_type_fasta.txt'
    PHOSPHATE_LIB_FASTA_FILE_ftp_link <- stringr::str_replace_all(PHOSPHATE_LIB_FASTA_FILE_ftp_link, 'fasta_type', fasta_type)
    PHOSPHATE_LIB_FASTA_FILE_ftp_link <- stringr::str_replace_all(PHOSPHATE_LIB_FASTA_FILE_ftp_link, 'species', species)
    PHOSPHATE_LIB_FASTA_FILE_data_type <- 'txt'
    PHOSPHATE_LIB_FASTA_DATA <- load_data_with_ftp(PHOSPHATE_LIB_FASTA_FILE_ftp_link, PHOSPHATE_LIB_FASTA_FILE_data_type)
    message('Save ', fasta_type, ' fasta file of ', species, ' to ', PHOSPHATE_LIB_FASTA_FILE_PATH)
    write.table(PHOSPHATE_LIB_FASTA_DATA, PHOSPHATE_LIB_FASTA_FILE_PATH, row.names = FALSE, col.names = TRUE, sep = '\t')
    message('Save successfully.')
  }else{
    PHOSPHATE_LIB_FASTA_DATA <- utils::read.table(file=PHOSPHATE_LIB_FASTA_FILE_PATH, header=TRUE, sep="\t")
  }
  fasta_data <- PHOSPHATE_LIB_FASTA_DATA
  # colnames(fasta_data) <- c('ID', 'Sequence')
  ###########################################################################################################

  border_limit <- floor(fixed_length/2)
  aligned_seq <- NULL
  GI_nrow <- length(ID)
  cat('Pre-align:', GI_nrow, 'phos-pepitdes.\n')
  cat('Fixed sequence length is ', fixed_length, '.\n', sep = '')
  cat('It needs few time.\n')
  for(i in seq_len(GI_nrow)){
    gi <- ID[i]
    aa_index <- AA_in_protein[i]
    loc_index <- as.numeric(stringr::str_split(aa_index, "[STY]", n = Inf, simplify = FALSE)[[1]])[2]
    index <- which(fasta_data[,1] == gi)
    if(length(index) > 0){
      refseq <- as.vector(fasta_data[index,2])
      refseq_len <- nchar(refseq)

      left_limit <- loc_index - border_limit
      right_limit <- loc_index + border_limit

      if(left_limit>=1 & right_limit>refseq_len){
        right_limit <- refseq_len
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "right", pad = '_')
      }else if(left_limit<1 & right_limit<=refseq_len){
        left_limit <- 1
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "left", pad = '_')
      }else if(left_limit<1 & right_limit>refseq_len){
        left_limit <- 1
        right_limit <- refseq_len
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "both", pad = '_')
      }else{
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
      }
    }else{
      truncated_seq <- NA
    }
    aligned_seq <- c(aligned_seq, truncated_seq)
    if(i %% 5000 == 0){
      cat('Aligned:', i, 'phos-pepitdes.\n')
    }
    if(i == GI_nrow){
      cat('Aligned:', i, 'phos-pepitdes.\n')
      cat('Finish OK! ^_^\n')
    }

  }
  cat('\n')
  aligned_sequence_df_based_on_fasta_library <- data.frame(ID, Sequence, AA_in_protein, aligned_seq)
  return(aligned_sequence_df_based_on_fasta_library)
}
