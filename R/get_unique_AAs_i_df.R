#' Get data frame without redundancy.

#'
#' @param df_with_AAs_i a data frame for peptides of the ith protein.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame with sites in unique protein.
#' @export
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/get_unique_AAs_i_df.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "get_unique_AAs_i_df.RData")
#' load("get_unique_AAs_i_df.RData")
#'
#' summary_df_of_unique_protein_with_sites <- get_unique_AAs_i_df(df_with_AAs_i)
#' head(summary_df_of_unique_protein_with_sites)
#'
get_unique_AAs_i_df <- function(df_with_AAs_i){
  # 3
  # Remove redundant records for specific protein to keep psites in protein sequence unique.
  protein_sites <- as.vector(df_with_AAs_i$AA_in_protein)
  unique_protein_sites <- sort(unique(protein_sites))
  unique_protein_sites_count <- length(unique_protein_sites)
  summary_df_with_unique_protein_sites <- c()
  for(i in seq_len(unique_protein_sites_count)){
    unique_protein_site <- unique_protein_sites[i]
    index_of_mapping_unique_protein_site <- which(protein_sites==unique_protein_site)
    df_with_unique_protein_site <- df_with_AAs_i[index_of_mapping_unique_protein_site,]
    AA_in_protein <- unique_protein_site
    AA_in_peptide <- paste(df_with_unique_protein_site$AA_in_peptide, collapse = '||')
    Sequence <- paste(df_with_unique_protein_site$Sequence, collapse = '||')
    ID <- as.vector(df_with_unique_protein_site$ID[1])
    Modification <- as.vector(df_with_unique_protein_site$Modification[1])
    GeneSymbol <- paste(df_with_unique_protein_site$GeneSymbol, collapse = '||')
    Value_Sum <- colSums(df_with_unique_protein_site[,-c(seq_len(6))])
    summary_df_with_unique_protein_site <- cbind(data.frame(AA_in_protein, AA_in_peptide, Sequence, ID, Modification, GeneSymbol),
                                                 t(data.frame(Value_Sum)))
    summary_df_with_unique_protein_sites <- rbind(summary_df_with_unique_protein_sites, summary_df_with_unique_protein_site)
  }
  return(summary_df_with_unique_protein_sites)
}
