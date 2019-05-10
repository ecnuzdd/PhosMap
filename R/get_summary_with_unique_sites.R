
#' Assign psites to protein sequence.
#'
#' Construct the data frame with unique phosphorylation site for each protein sequence and eliminate redundancy.
#'
#'
#' @param combinated_df_with_mapped_gene_symbol A dataframe with Sequence, ID, Modification, Gene Symbol, Area and PSMs as input.
#' @param species A string, the options are human, mouse and rattus, the default is human.
#' @param fasta_type, A string for fasta source, the options are refseq and uniprot, the default is refseq
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe that all redundant psites are assigned to protein sequence.
#' @export
#'
#' @examples
#' \dontrun{
#' summary_df_of_unique_proteins_with_sites <- get_summary_with_unique_sites(
#'   combinated_df_with_mapped_gene_symbol,
#'   species = 'human',
#'   fasta_type = 'refseq'
#' )
#' }
get_summary_with_unique_sites <- function(
  combinated_df_with_mapped_gene_symbol,
  species = 'human',
  fasta_type = 'refseq'
){
  requireNamespace('utils')
  requireNamespace('stringr')
  # unique phosphorylation sites
  cat('\n The 6th step: construct the data frame with unique phosphorylation site for each protein sequence.')

  ###########################################################################################################
  # Read fasta file: fasta_data
  # fasta_type <- 'refseq' # 'uniprot'
  if(fasta_type == 'refseq' | fasta_type == 'uniprot'){
    PHOSPHATE_LIB_FASTA_DIR <- normalizePath(
      system.file(
        'extdata', 'fasta_libarary',
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

  # https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/fasta_libarary/refseq/human/human_ref_fasta.txt
  if(!file.exists(PHOSPHATE_LIB_FASTA_FILE_PATH)){
    PHOSPHATE_LIB_FASTA_FILE_http_link <- 'https://raw.githubusercontent.com/ecnuzdd/PhosMap_datasets/master/fasta_libarary/fasta_type/species/species_fasta_type_fasta.txt'
    PHOSPHATE_LIB_FASTA_FILE_http_link <- stringr::str_replace_all(PHOSPHATE_LIB_FASTA_FILE_http_link, 'fasta_type', fasta_type)
    PHOSPHATE_LIB_FASTA_FILE_http_link <- stringr::str_replace_all(PHOSPHATE_LIB_FASTA_FILE_http_link, 'species', species)
    PHOSPHATE_LIB_FASTA_FILE_data_type <- 'txt'
    PHOSPHATE_LIB_FASTA_DATA <- load_data_with_http(PHOSPHATE_LIB_FASTA_FILE_http_link, PHOSPHATE_LIB_FASTA_FILE_data_type)
    message('Save ', fasta_type, ' fasta file of ', species, ' to ', PHOSPHATE_LIB_FASTA_FILE_PATH)
    write.table(PHOSPHATE_LIB_FASTA_DATA, PHOSPHATE_LIB_FASTA_FILE_PATH, row.names = FALSE, col.names = TRUE, sep = '\t')
    message('Save successfully.')
  }else{
    PHOSPHATE_LIB_FASTA_DATA = utils::read.table(file=PHOSPHATE_LIB_FASTA_FILE_PATH, header=TRUE, sep="\t")
  }
  fasta_data <- PHOSPHATE_LIB_FASTA_DATA
  # colnames(fasta_data) <- c('ID', 'Sequence')
  ###########################################################################################################




  id_data <- combinated_df_with_mapped_gene_symbol

  # Keep peptides assigned to unique protein
  id_data_only_peptide2gi <- id_data[which(!grepl(';', as.vector(id_data$ID))),]

  # Determine locations of the psites each peptide mapped to protein squence.
  modification_index_in_protein_seq_list <- get_modification_index(id_data_only_peptide2gi,
                                                                   fasta_data)

  proteins_in_id_data_only_peptide2gi <- as.vector(id_data_only_peptide2gi$ID)
  sequences_in_id_data_only_peptide2gi <- as.vector(id_data_only_peptide2gi$Sequence)
  value_in_id_data_only_peptide2gi <- id_data_only_peptide2gi[, -c(seq_len(4))]

  unique_proteins <- unique(proteins_in_id_data_only_peptide2gi)
  unique_protein_count <- length(unique_proteins)

  # Show psites and modifications of one protein, merge the values with the same modification type.
  cat('\n', 'Map phosphorylation sites to protein sequence and eliminate redundancy.')
  system.time({
    summary_df_of_unique_proteins_with_sites <- c()
    for(i in seq_len(unique_protein_count)){

      df_with_AAs_i <- get_df_with_AAs_i(unique_proteins,
                                         i,
                                         id_data_only_peptide2gi,
                                         proteins_in_id_data_only_peptide2gi,
                                         sequences_in_id_data_only_peptide2gi,
                                         modification_index_in_protein_seq_list)


      summary_df_of_unique_protein_with_sites <- get_unique_AAs_i_df(df_with_AAs_i)


      summary_df_of_unique_proteins_with_sites <- rbind(
        summary_df_of_unique_proteins_with_sites,
        summary_df_of_unique_protein_with_sites
      )

      if(i%%500==0 | i == unique_protein_count){
        cat('\n completed: ', i, '/', unique_protein_count)
      }
    }
  })
  summary_df_of_unique_proteins_with_sites_rownames <- paste(as.vector(summary_df_of_unique_proteins_with_sites$ID),
                                                             as.vector(summary_df_of_unique_proteins_with_sites$AA_in_protein),
                                                             sep = '_')
  rownames(summary_df_of_unique_proteins_with_sites) <- summary_df_of_unique_proteins_with_sites_rownames
  summary_df_of_unique_proteins_with_sites_colnames <- colnames(summary_df_of_unique_proteins_with_sites)
  index_of_PSMs <- which(grepl('_PSMs', summary_df_of_unique_proteins_with_sites_colnames))
  if(length(index_of_PSMs)>0){
    summary_df_of_unique_proteins_with_sites <- summary_df_of_unique_proteins_with_sites[,-index_of_PSMs]
  }
  summary_df_of_unique_proteins_with_sites$GeneSymbol <- apply(data.frame(summary_df_of_unique_proteins_with_sites$GeneSymbol),
                                                               1,
                                                               function(x){
                                                                 if(grepl('||', x)){
                                                                   x <- as.vector(x)
                                                                   x <- strsplit(x, split = '||', fixed = TRUE)
                                                                   x[[1]][1]
                                                                 }
                                                               })
  cat('\n The 6th step: construct over.')
  return(summary_df_of_unique_proteins_with_sites)
}
