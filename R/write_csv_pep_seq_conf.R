#' Write data to specific direction with CSV format.
#'
#' @param expName a string for experiment name as input.
#' @param outputName a string for experiment name as output.
#' @param mascotfileNames a vector for storing mascot file names.
#' @param refFileName a string for reference file name.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @return Write data to specific direction with CSV format.
#' @export
#' @examples
#' \dontrun{
#' write_csv_pep_seq_conf(expName,
#'   outputName, mascot_txt_dir_path_expName_path,
#'   firmiana_peptide_dir_path_expName_path)
#' }

write_csv_pep_seq_conf <- function(expName, outputName, mascotfileNames, refFileName){
  requireNamespace('utils')
  fileNames_len <- length(mascotfileNames)
  fileData_list <- list()
  for(i in seq_len(fileNames_len)){
    mascotfileName <- mascotfileNames[i]
    df <- get_filtered_df(mascotfileName, refFileName)
    fileData_list[[i]] <- df
  }
  mergeMat <- subset(fileData_list[[1]], select = c('pep_seq', 'pep_var_mod_conf'))
  if(fileNames_len>1){
    for(i in 2:fileNames_len){
      tmp_mergeMat <- subset(fileData_list[[i]], select = c('pep_seq', 'pep_var_mod_conf'))
      mergeMat <- merge(mergeMat, tmp_mergeMat, by='pep_seq', all = TRUE)
    }
  }
  mergeMat_ionscore <- subset(fileData_list[[1]], select = c('pep_seq', 'pep_score'))
  if(fileNames_len>1){
    for(i in 2:fileNames_len){
      tmp_mergeMat_ionscore <- subset(fileData_list[[i]], select = c('pep_seq', 'pep_score'))
      mergeMat_ionscore <- merge(mergeMat_ionscore, tmp_mergeMat_ionscore, by='pep_seq', all = TRUE)
    }
  }
  mergeMat_row <- nrow(mergeMat)
  pep_seq <- as.vector(mergeMat$pep_seq)
  pep_var_mod_conf <- NULL
  for(i in seq_len(mergeMat_row)){
    score <- as.vector(unlist(mergeMat_ionscore[i,-1]))
    conf <- as.vector(unlist(mergeMat[i,-1]))
    index <- which.max(score)
    conf <- conf[index]
    pep_var_mod_conf <- c(pep_var_mod_conf, conf)
  }
  df1 <- data.frame(pep_seq, pep_var_mod_conf)
  utils::write.csv(df1, outputName, row.names = FALSE)
}
