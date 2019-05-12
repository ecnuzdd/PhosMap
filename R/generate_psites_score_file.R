#' Generate peptide identification files with psites scores.
#'
#' Based on mascot txt files with psites and peptide identification files downloaded from Firmiana, the file with phosphorylation modifications is generated.
#'
#'
#' @param mascot_txt_dir A folder containing identification xml files with psites scores as input.
#' @param firmiana_peptide_dir A folder containing identification txt files downloaded from Firmiana as input.
#' @param psites_score_dir A folder used for saving files of peptide identification files with psites scores
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A series of output files saved in the psites_score_dir
#' @export
#'
#' @examples
#' \dontrun{
#' generate_psites_score_file(mascot_txt_dir, firmiana_peptide_dir, psites_score_dir)
#' }
#'

generate_psites_score_file <- function(mascot_txt_dir, firmiana_peptide_dir, psites_score_dir){
  mascot_txt_dir_paths <- list.dirs(mascot_txt_dir)
  mascot_txt_dir_paths <- mascot_txt_dir_paths[-1]
  mascot_txt_dir_paths_expNames <- list.files(mascot_txt_dir)

  firmiana_txt_file_names <- list.files(firmiana_peptide_dir)
  firmiana_txt_file_names_expNames <- apply(data.frame(firmiana_txt_file_names), 1, function(x){
    x <- strsplit(x, split = '_')[[1]]
    x[1]
  })

  mascot_txt_dir_paths_len <- length(mascot_txt_dir_paths)
  cat('\n Total file: ', mascot_txt_dir_paths_len)
  cat('\n It will take a little while.')
  for(i in seq_len(mascot_txt_dir_paths_len)){
    mascot_txt_dir_path <- mascot_txt_dir_paths[i]
    mascot_txt_dir_path_expName <- mascot_txt_dir_paths_expNames[i]
    mascot_txt_dir_path_expName_path <- list.files(mascot_txt_dir_path)
    mascot_txt_dir_path_expName_path <- normalizePath(
      file.path(mascot_txt_dir_path, mascot_txt_dir_path_expName_path)
    )

    match_index <- match(mascot_txt_dir_path_expName, firmiana_txt_file_names_expNames)
    firmiana_txt_file_name <- firmiana_txt_file_names[match_index]
    firmiana_peptide_dir_path_expName_path <- normalizePath(
      file.path(firmiana_peptide_dir, firmiana_txt_file_name)
    )

    outputName <- normalizePath(
      file.path(psites_score_dir, paste(mascot_txt_dir_path_expName, '_psites_score.csv', sep = '')),
      mustWork = FALSE
    )
    expName <- mascot_txt_dir_path_expName
    write_csv_pep_seq_conf(expName, outputName, mascot_txt_dir_path_expName_path, firmiana_peptide_dir_path_expName_path)

    cat('\n Completed file: ', i, '/', mascot_txt_dir_paths_len)

  }

}
