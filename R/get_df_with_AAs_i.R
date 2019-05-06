#' Get data frame of amino acide sequences for a protein.

#'
#' @param unique_proteins a vector for unique proteins.
#' @param i the ith unique proteins.
#' @param id_data_only_peptide2gi a data frame for peptides with protein gi.
#' @param proteins_in_id_data_only_peptide2gi a vector for proteins with only protein gi.
#' @param sequences_in_id_data_only_peptide2gi a vector for peptides with only protein gi.
#' @param modification_index_in_protein_seq_list a list for the index of modifications in protein sequence.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame with sequences for a protein.
#' @export
#'


get_df_with_AAs_i <- function(unique_proteins,
                              i,
                              id_data_only_peptide2gi,
                              proteins_in_id_data_only_peptide2gi,
                              sequences_in_id_data_only_peptide2gi,
                              modification_index_in_protein_seq_list) {
  # 2
  # Get the phosphorylation sequence matrix from the ist unique protein
  unique_protein <- unique_proteins[i]
  index_of_mapping <- which(proteins_in_id_data_only_peptide2gi == unique_protein) # All peptides location (in id_data_only_peptide2gi) assigned to one protein.
  index_of_mapping_count <- length(index_of_mapping)

  sequences_of_mapping <- sequences_in_id_data_only_peptide2gi[index_of_mapping] # All peptide sequence assigned to one protein.
  df_of_mapping <- id_data_only_peptide2gi[index_of_mapping, ] # All peptide sequence and related information assigned to one protein.

  df_with_AAs_i <- c()
  for (j in seq_len(index_of_mapping_count)) {
    index_of_mapping_j <- index_of_mapping[j]

    # Store psites location index of all peptides assigned to one protein which are mapped to the assigned protein sequence.
    site_index_in_protein_seq <- modification_index_in_protein_seq_list[[index_of_mapping_j]]

    # Store psites location index of all peptides assigned to one protein which are mapped to the assigned peptide sequence.
    sequences_of_mapping_j <-  sequences_of_mapping[j]
    site_index_in_peptide_seq <- unlist(gregexpr("[a-z]", sequences_of_mapping_j))

    df_of_mapping_j <- df_of_mapping[j, ]

    modified_AA_count_j <- length(site_index_in_peptide_seq)
    AAs_in_peptide <- NULL # Store psites and location index of all peptides assigned to one protein which are mapped to the assigned protein sequence.
    AAs_in_protein <- NULL # Store psites and location index of all peptides assigned to one protein which are mapped to the assigned peptide sequence.
    df_with_AAs_j <- c()

    for (k in seq_len(modified_AA_count_j)) {
      site_index_in_peptide_seq_k <- site_index_in_peptide_seq[k]
      site_index_in_protein_seq_k <- site_index_in_protein_seq[k]
      AA <- substr(
        sequences_of_mapping_j,
        site_index_in_peptide_seq_k,
        site_index_in_peptide_seq_k
      )
      AA_in_peptide <- paste(AA, site_index_in_peptide_seq_k, sep = '')
      AA_in_protein <- paste(AA, site_index_in_protein_seq_k, sep = '')
      AAs_in_peptide <- c(AAs_in_peptide, AA_in_peptide)
      AAs_in_protein <- c(AAs_in_protein, AA_in_protein)
      df_with_AAs_j_i <- data.frame(AA_in_protein, AA_in_peptide, df_of_mapping_j)
      df_with_AAs_j <- rbind(df_with_AAs_j, df_with_AAs_j_i)
    }
    df_with_AAs_i <- rbind(df_with_AAs_i, df_with_AAs_j)


  }
  return(df_with_AAs_i)
}
