#' Get peptides data frame passed phosphorylation sites quality control.
#'
#' Filter phosphorylation sites by extracting all peptides with ion score>=20 and FDR<0.01 from Firmiana and having psites scores.
#' Generate new IDs consisting of sequence, gi, psite.
#' Quantification values containing area and psm.
#'
#' @param firmiana_peptide_dir A folder containing peptide identification files from Firmiana as input.
#' @param psites_score_dir A folder containing psites scores files extracted from mascot xml as input.
#' @param phospho_experiment_design_file_path A string representing the path of phospho-experiment design file as input.
#' @param qc A boolean value representing whether it has QC files. The default is True.
#' @param min_score A numeric for the minimum score of credible peptides, the default is 20 for Mascot ion score.
#' @param min_FDR A numeric for the minimum FDR of credible peptides, the default is 0.01.
#'
#' @return A merged data frame containing sequence, gi, psite, area and psm.
#' @export
#'
#' @examples
#' \dontrun{
#' merge_df_with_phospho_peptides = pre_process_filter_psites(
#'   firmiana_peptide_dir,
#'   psites_score_dir
#' )
#' }
pre_process_filter_psites <- function(firmiana_peptide_dir, psites_score_dir,
                                      phospho_experiment_design_file_path, qc,
                                      min_score = 20, min_FDR = 0.01) {
  requireNamespace('utils')
  PEPTIDE_DIR <- normalizePath(firmiana_peptide_dir, mustWork = FALSE)
  if(!file.exists(firmiana_peptide_dir)){
    cat(firmiana_peptide_dir, ' -> ', 'No the directory.')
    stop('')
  }

  PSITES_WITH_SCORE_DIR <- normalizePath(psites_score_dir, mustWork = FALSE)
  if(!file.exists(psites_score_dir)){
    cat(psites_score_dir, ' -> ', 'No the directory.')
    stop('')
  }

  phospho_experiment_design_file_path <- normalizePath(phospho_experiment_design_file_path, mustWork = FALSE)
  if(!file.exists(phospho_experiment_design_file_path)){
    cat(phospho_experiment_design_file_path, ' -> ', 'No the file')
    stop('')
  }

  # read experiment design file and make merged experments keep order of experiment design
  phospho_experiment_design_file <- utils::read.table(phospho_experiment_design_file_path, sep = '\t',
                                               header = TRUE, stringsAsFactors = NA)
  phospho_experiment_ID <- as.vector(unlist(phospho_experiment_design_file$Experiment_Code))

  #### (1) Read peptide identification file from Firmiana ####
  cat('\n The 1st step: read peptide identification files.')
  result_list_from_PEPTIDE_DIR <- get_file_info_from_dir(PEPTIDE_DIR, phospho_experiment_ID)
  files <- result_list_from_PEPTIDE_DIR$file_data_list
  peptide.id <- result_list_from_PEPTIDE_DIR$file_ID


  #### (2) Read psits score file ####
  if(qc){
    cat('\n The 2nd step: read psites QC files.')
    result_list_from_PSITES_WITH_SCORE_DIR <- get_file_info_from_dir(PSITES_WITH_SCORE_DIR,
                                                                                          phospho_experiment_ID)
    files_site_score <- result_list_from_PSITES_WITH_SCORE_DIR$file_data_list
    site_score.id <- result_list_from_PSITES_WITH_SCORE_DIR$file_ID
  }else{
    cat('\n The 2nd step: no QC files.')
  }



  #### (3) Filter sites ####
  # Filter phosphorylation sites by extracting all peptide identification with ion score>=20 and FDR<0.01 from Firmiana and having psites score.
  # ************
  # *required column:
  # *file_peptide: Ion.Score, FDR, Area, PSMs, Sequence, Protein.Groups.Accessions, Modification
  # *file_site_score: pep_seq, pep_var_mod_conf
  cat('\n The 3rd step: filter peptides based on site quality.')
  if(qc){
    result_list_with_filtered_sites <- get_list_with_filted_sites(peptide.id, files,
                                                                 files_site_score, qc,
                                                                 min_score, min_FDR)
  }else{
    result_list_with_filtered_sites <- get_list_with_filted_sites(peptide.id, files,
                                                                 NULL, qc,
                                                                 min_score, min_FDR)
  }

  peptide_df_with_area_psm_list <- result_list_with_filtered_sites$peptide_df_with_area_psm_list # including: area, psm
  ID_of_seq_gi_site_list <- result_list_with_filtered_sites$ID_of_seq_gi_site_list # including: seq_gi_psite
  ID_DF_list <- result_list_with_filtered_sites$ID_DF_list # including: seq_gi_psite, area, psm


  #### (4) Based on unique peptide, merge all experiments ####
  cat('\n The 4th step: merge data based on peptides (unique ID).')
  merge_df_with_phospho_peptides <- get_merged_phospho_df(peptide.id,
                                                          peptide_df_with_area_psm_list,
                                                          ID_of_seq_gi_site_list, ID_DF_list)

  # delete psm column
  merge_df_with_phospho_peptides_colnames <- colnames(merge_df_with_phospho_peptides)
  index_of_PSMs <- grep('_PSMs', merge_df_with_phospho_peptides_colnames)
  merge_df_with_phospho_peptides <- merge_df_with_phospho_peptides[,-index_of_PSMs]



  merge_df_with_phospho_peptides_colnames <- colnames(merge_df_with_phospho_peptides)
  ID <- as.vector(merge_df_with_phospho_peptides[,1])
  Value <- merge_df_with_phospho_peptides[,-1]
  Value_colnames <- colnames(Value)
  Value_colnames_ID <- apply(data.frame(Value_colnames), 1, function(x){
    x <- strsplit(x, split = '_')[[1]][1]
    x
  })
  index_of_match <- match(phospho_experiment_ID, Value_colnames_ID)
  Value <- Value[,index_of_match]
  merge_df_with_phospho_peptides <- data.frame(ID, Value)
  colnames(merge_df_with_phospho_peptides) <- c(merge_df_with_phospho_peptides_colnames[1], phospho_experiment_ID)
  return(merge_df_with_phospho_peptides)
}







