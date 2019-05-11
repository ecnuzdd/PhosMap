#' Get merged data frame with phospho-peptides.
#'
#' @param peptide_id a vector for peptide ID.
#' @param peptide_df_with_area_psm_list a list for peptides with areas and PSMs.
#' @param ID_of_seq_gi_site_list a list for peptides ID with sequence, gi and site.
#' @param ID_DF_list a list for ID and values.
#'
#' @return a merged data frame with phospho-peptides.
#' @export

get_merged_phospho_df <- function(
  peptide_id,
  peptide_df_with_area_psm_list,
  ID_of_seq_gi_site_list,
  ID_DF_list
){
  cat('\n The 4th step is running.')
  # IDs, accession, modification identified in all files.
  peptides_in_all_experiments <- unique(unlist(ID_of_seq_gi_site_list))
  # Reserve peptides only phosphorylation modified.
  phospho_peptides_in_all_experiments <- peptides_in_all_experiments[grepl('Phospho', peptides_in_all_experiments)]
  peptide_id_len = length(peptide_id)
  df_list_with_phospho_peptides <- list()
  for(i in seq_len(peptide_id_len)){
    ID_of_seq_gi_site <- ID_of_seq_gi_site_list[[i]]
    peptide_df_with_area_psm <- peptide_df_with_area_psm_list[[i]]
    index_of_match_phospho_peptides <- match(ID_of_seq_gi_site, phospho_peptides_in_all_experiments)
    index_of_NonNA <- which(!is.na(index_of_match_phospho_peptides))
    df_with_phospho_peptides <- ID_DF_list[[i]][index_of_NonNA,]
    df_list_with_phospho_peptides[[i]] <- df_with_phospho_peptides
  }
  merge_df_with_phospho_peptides <- df_list_with_phospho_peptides[[1]]
  if(peptide_id_len>1){
    for(i in 2:peptide_id_len){
      tmp_df <- df_list_with_phospho_peptides[[i]]
      merge_df_with_phospho_peptides <- merge(merge_df_with_phospho_peptides, tmp_df, by = 'ID', all = TRUE)
    }
  }
  ID_of_seq_gi_site <- as.vector(merge_df_with_phospho_peptides$ID)
  merge_df_with_phospho_peptides_Value <- as.matrix(merge_df_with_phospho_peptides[,-1])
  index_of_NA <- which(is.na(merge_df_with_phospho_peptides_Value))
  if(length(index_of_NA)>0){
    merge_df_with_phospho_peptides_Value[index_of_NA] <- 0
  }
  merge_df_with_phospho_peptides <- data.frame(ID_of_seq_gi_site, merge_df_with_phospho_peptides_Value)
  cat('\n The 4th step is over ^_^.')
  return(merge_df_with_phospho_peptides)
}
