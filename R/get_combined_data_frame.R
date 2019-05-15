#' Get a data frame mapped ID to Gene Symbol.
#'
#' This is an intermediate file and a dataframe with Gene Symbol exported.
#' Based on a library file consisting of mapping relationships about Gene Symbol, GeneID, RefSeq_Protein_GI, RefSeq_Protein_Accession and Uniprot_Protein_Accession,
#' a new dataframe with Sequence, GI, Modification, Gene Symbol, Area and PSMs,is contructed.
#'
#' @param merge_df_with_phospho_peptides A dataframe consisting of IDs (Sequence_GI_Psite) and Area values.
#' @param species A string, the options are human, mouse and rattus, the default is human.
#' @param id_type A string, the options are 'GeneID', 'RefSeq_Protein_GI', 'RefSeq_Protein_Accession' and 'Uniprot_Protein_Accession', the default is RefSeq_Protein_GI.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe with Sequence, GI, Modification, Gene Symbol, Area values and PSMs
#' @export
#'
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/merge_df_with_phospho_peptides.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "merge_df_with_phospho_peptides.RData")
#' load("merge_df_with_phospho_peptides.RData")
#'
#' combined_df_with_mapped_gene_symbol <- get_combined_data_frame(
#'   merge_df_with_phospho_peptides, species = 'human',
#'   id_type = 'RefSeq_Protein_GI')
#'
#' head(combined_df_with_mapped_gene_symbol)
#'

get_combined_data_frame <- function(
  merge_df_with_phospho_peptides,
  species = 'human',
  id_type = 'RefSeq_Protein_GI'
){
  # Read library file, map GI to Gene Symbol
  requireNamespace('utils')
  requireNamespace('stringr')

  cat('\n The 5th step: write the data frame with symbols mapping to genes.')

  ######################################################################################
  # load datasets
  id_coversion_table_dir <- normalizePath(
    system.file(
      'extdata',
      'id_coversion_table',
      package = "PhosMap"
    ),
    mustWork = FALSE
  )

  PHOSPHATE_LIB_MAPPING_FILE_PATH <- normalizePath(
    file.path(id_coversion_table_dir, paste(species, 'ID.txt', sep = '_')),
    mustWork = FALSE
  )

  if(!file.exists(PHOSPHATE_LIB_MAPPING_FILE_PATH)){
    id_coversion_table_ftp_link <- 'ftp://111.198.139.72:4000/pub/PhosMap_datasets/id_coversion_table/species_ID.txt'
    id_coversion_table_ftp_link <- stringr::str_replace_all(id_coversion_table_ftp_link, 'species', species)
    id_coversion_table_data_type <- 'txt'
    id_coversion_table <- load_data_with_ftp(id_coversion_table_ftp_link, id_coversion_table_data_type)
    message('Save id coversion table of ', species, ' to ', PHOSPHATE_LIB_MAPPING_FILE_PATH)
    # write.csv(id_coversion_table, PHOSPHATE_LIB_MAPPING_FILE_PATH, row.names = FALSE)
    write.table(id_coversion_table, PHOSPHATE_LIB_MAPPING_FILE_PATH, sep = '\t', row.names = FALSE)
    message('Save successfully.')
  }else{
    id_coversion_table = utils::read.table(PHOSPHATE_LIB_MAPPING_FILE_PATH, sep = '\t', header = TRUE)
  }
  ######################################################################################


  cat('\n The 5th step is running.')
  # Split a string: sequenceID, accession, modification
  seq_gi_site_vector <- as.vector(merge_df_with_phospho_peptides$ID_of_seq_gi_site)
  Sequence <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][1]
  })
  ID <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][2]
  })
  Modification <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][3]
  })


  ##########################################################################################################
  # id_types <- c('GeneID', 'RefSeq_Protein_GI', 'RefSeq_Protein_Accession', 'Uniprot_Protein_Accession')
  # GeneSymbol
  # construct dict
  id_type <- 'RefSeq_Protein_GI'
  MappingDf <- id_coversion_table[, c('GeneSymbol', id_type)]
  invalid_index <- which(as.vector(unlist(MappingDf[,2])) == '' | as.vector(unlist(MappingDf[,2])) == '-')
  if(length(invalid_index)>0){
    MappingDf <- MappingDf[-invalid_index,]
  }
  MappingDf_row <- nrow(MappingDf)
  cat('\n', 'Construct dictionary based on GeneSymbol and specific ID.')
  mapping_dict <- NULL
  cat('\n', 'The total:', MappingDf_row)
  for(i in seq_len(MappingDf_row)){
    x <- as.vector(MappingDf[i,1])
    y <- as.vector(unlist(MappingDf[i,2]))
    y <- strsplit(y, split = '; ')[[1]]
    x_v <- rep(x, length(y))
    names(x_v) <- y
    mapping_dict <- c(mapping_dict, x_v)
    if(i%%5000==0 | i == MappingDf_row){
      cat('\n', 'Completed:', i, '/', MappingDf_row)
    }
  }
  ##########################################################################################################

  GeneSymbol <- apply(data.frame(ID), 1, function(x, mapping_dict, id_type){
    gi_all <- strsplit(x, split=";", fixed = TRUE)[[1]]

    gi_mapping_symbol <- apply(data.frame(gi_all), 1, function(y, mapping_dict, id_type){
      if(id_type == 'RefSeq_Protein_GI'){
        y = stringr::str_replace_all(y, 'gi[|]', '')
      }
      return(mapping_dict[y])
    }, mapping_dict = mapping_dict, id_type)

    gi_mapping_symbol_unique <- unique(gi_mapping_symbol[which(!is.na(gi_mapping_symbol))])
    gi_mapping_symbol_unique_count <- length(gi_mapping_symbol_unique)


    if(gi_mapping_symbol_unique_count == 0){
      return(NA)
    }else if(gi_mapping_symbol_unique_count == 1){
      return(gi_mapping_symbol_unique)
    }else{
      return(paste(gi_all, collapse = ';'))
    }
  }, mapping_dict = mapping_dict, id_type = id_type)


  # sequenceID, accession, symbol, modification, quantification_value_in_experiment
  df_of_combination <- data.frame(Sequence, ID, Modification, GeneSymbol, merge_df_with_phospho_peptides[,-1]) # delete first column
  index_of_NonNA <- which(!is.na(GeneSymbol))
  df_of_combination <- df_of_combination[index_of_NonNA,]
  cat('\n The 5th step is over ^_^.')
  cat('\n The 5th step: write the data frame with symbols mapping to genes.')
  return(df_of_combination)
}
