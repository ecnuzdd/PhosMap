#' Filter phosphorylation sites.
#'
#' Filter phosphorylation sites by extracting all peptides with ion score>=20 and FDR<0.01 from Firmiana and having psites scores.
#'
#' @param peptide_id A vector containing experiment ids as input.
#' @param files A data list containing peptides identificaton from Firmiana as input.
#' @param files_site_score A data list containing psites scores extracted from mascot xml. The default is NULL,
#' which represents no QC file.
#' @param qc A boolean value representing whether it has QC files. The default is True.
#' @param min_score A numeric for the minimum score of credible peptides, the default is 20 for Mascot ion score.
#' @param min_FDR A numeric for the minimum FDR of credible peptides, the default is 0.01.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A list containing peptides dataframe with area values and psm,
#' IDs with mergered sequences, gi and sites, new peptides dataframe combined previous peptides dataframe and IDs.
#' @export
#'
#' @examples
#' \dontrun{
#' result_list_with_filtered_sites = get_list_with_filtered_sites(
#'   peptide.id,
#'   files,
#'   files_site_score
#' )
#' }
get_list_with_filtered_sites <- function(
  peptide_id,
  files,
  files_site_score,
  qc,
  min_score,
  min_FDR
){
  peptide_df_with_area_psm_list <- list() # data.frame(area, psm)
  ID_of_seq_gi_site_list <- list() # seq_gi_psite
  ID_DF_list <- list() # seq_gi_psite + data.frame(area, psm)
  peptide_id_len <- length(peptide_id) # File Numbers
  # ************
  # *Required column:
  # *file_peptide: Ion_Score, FDR, Area, PSMs, Sequence, Protein_Groups_Accessions, Modification
  # *file_site_score: pep_seq, pep_var_mod_conf
  cat('\n Total file: ', peptide_id_len)
  for(i in seq_len(peptide_id_len)){
    cat('\n completed: ',i,'/',peptide_id_len)

    file_peptide <- data.frame(files[[i]])
    # Set parameters 1：reserve peptides with ion score><-20 and FDR<0.01.
    index_of_row_filters_meet_ionscore_and_FDR <- which(as.numeric(as.vector(file_peptide$Ion.Score)) >= min_score &
                                                         as.numeric(as.vector(file_peptide$FDR)) < min_FDR)
    file_peptide <- file_peptide[index_of_row_filters_meet_ionscore_and_FDR, ]

    if(!qc){
      file_peptide_subset <- file_peptide
    }else{
      # Extract peptides with psites score.
      file_site_score  <-  as.data.frame(files_site_score[[i]])
      index_of_row_filters_have_site_score <- which(grepl('%', file_site_score$pep_var_mod_conf))
      file_site_score  <-  file_site_score[index_of_row_filters_have_site_score,]

      # Reserve peptides with psites score in file_peptide.
      index_of_peptide_with_site_score_in_file_peptide <- match(as.vector(file_site_score[,1]), as.vector(file_peptide[,1]))
      index_of_NA <- which(is.na(index_of_peptide_with_site_score_in_file_peptide))
      if(length(index_of_NA)>0){
        index_of_peptide_with_site_score_in_file_peptide <- index_of_peptide_with_site_score_in_file_peptide[-index_of_NA]
      }
      file_peptide_subset <- file_peptide[index_of_peptide_with_site_score_in_file_peptide,]
    }
    area <- as.numeric(as.vector(file_peptide_subset$Area))
    psms <- as.numeric(as.vector(file_peptide_subset$PSMs))

    peptide_df_with_area_psm <- data.frame(area, psms)
    peptide_df_with_area_psm_colnames <- paste(peptide_id[i], c('Area', 'PSMs'), sep = '_')
    colnames(peptide_df_with_area_psm) <- peptide_df_with_area_psm_colnames

    sequence_id <- as.vector(file_peptide_subset$Sequence)
    accession <- as.vector(file_peptide_subset$Protein.Groups.Accessions)
    modification <- as.vector(file_peptide_subset$Modification)
    ID_of_seq_gi_site <- paste(sequence_id, accession, modification, sep = '||')

    ID_DF <- data.frame(ID_of_seq_gi_site, peptide_df_with_area_psm)
    colnames(ID_DF) <- c("ID", peptide_df_with_area_psm_colnames)


    peptide_df_with_area_psm_list[[i]] <- peptide_df_with_area_psm # area, psm
    ID_of_seq_gi_site_list[[i]] <- ID_of_seq_gi_site # seq_gi_psite
    ID_DF_list[[i]] <- ID_DF # seq_gi_psite, area, psm

  }
  result_list <- list(
    peptide_df_with_area_psm_list = peptide_df_with_area_psm_list,
    ID_of_seq_gi_site_list = ID_of_seq_gi_site_list,
    ID_DF_list = ID_DF_list
  )
  return(result_list)
}
