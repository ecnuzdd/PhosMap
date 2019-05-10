#' Create R code to call python for parsering mascot xml.
#'
#' Extract the confidence probability of phosphorylation sites (psites) from mascot xml.
#' One file containing experiment codes and one folder containing mascot xml as input, the another folder is required as output.
#' Python is required and the corresponding xml package is also required.
#'
#'
#' @param phosphorylation_exp_design_info_file_path A string representing the file path of experiment code, for examples: experiment_code.txt
#' @param mascot_xml_dir A folder containing identification xml files searched by Mascot as input, for examples: Exp020901_F1_R1.xml
#' @param mascot_txt_dir A folder used for saving files which contains the confidence of phosphorylation sites, for examples: Exp020901_F1_R1.txt
#'
#' @return A series of output file saved in the mascot_txt_dir
#' @export
#'
#' @examples
#' \dontrun{
#' extract_psites_score(
#'   phosphorylation_exp_design_info_file_path,
#'   mascot_xml_dir,
#'   mascot_txt_dir
#' )
#' }

extract_psites_score <- function(
  phosphorylation_exp_design_info_file_path,
  mascot_xml_dir,
  mascot_txt_dir
){
  requireNamespace('utils')
  phosphorylation_exp_design_info_file_path <- normalizePath(phosphorylation_exp_design_info_file_path)
  if (!file.exists(phosphorylation_exp_design_info_file_path)) {
    cat('\n', phosphorylation_exp_design_info_file_path, ' -> ', 'No the file.')
    stop('')
  }
  mascot_xml_dir <- normalizePath(mascot_xml_dir)
  if (!file.exists(mascot_xml_dir)) {
    cat('\n', mascot_xml_dir, ' -> ', 'No the directory.')
    stop('')
  }
  mascot_xml_dir_files <- list.files(mascot_xml_dir)

  mascot_txt_dir <- normalizePath(mascot_txt_dir)
  if (!file.exists(mascot_txt_dir)) {
    cat('\n', mascot_txt_dir, ' -> ', 'No the directory, create it.')
    dir.create(mascot_txt_dir)
  }

  command <- "python"
  path2script <- system.file("src", "XMLParser_mascot_dat.py", package = "PhosMap") # The location of python script called

  # path2script <- "w:/R/R-3.3.2/library/PhosMap/src/XMLParser_mascot_dat.py"
  path2script <- normalizePath(path2script, mustWork = FALSE)

  # Get experiments codes by reading txt files
  experiment_code <- utils::read.table(phosphorylation_exp_design_info_file_path,
                                       sep = '\t',
                                       header = TRUE)
  experiment_code <- as.vector(unlist(experiment_code$Experiment_Code))

  # match txt files to mascot_xml_dir
  experiment_match_index <- match(experiment_code, mascot_xml_dir_files)
  na_index <- which(is.na(experiment_match_index))
  if(length(na_index)>0){
    na_experiments <- experiment_code[na_index]
    cat('\n', 'The following experiments do not exist in', mascot_xml_dir, '\n')
    for(na_experiment in na_experiments){
      cat('\n', na_experiment, '\n')
    }
    stop('')
  }

  experiment_code_count <- length(experiment_code)
  if (experiment_code_count < 1) {
    cat('\n', phosphorylation_exp_design_info_file_path, '\n')
    stopifnot('No experiments')
  }

  cat('\n Start extracting the confidence of Psites from mascot.xml.')
  cat('\n Total ', experiment_code_count, ' experiment(s).')
  cat('\n It will take a little while.')

  parent_dir <- dirname(phosphorylation_exp_design_info_file_path)
  parent_dir <- normalizePath(parent_dir)
  log_dir <- normalizePath(file.path(parent_dir, 'log'), mustWork = FALSE)
  if (!file.exists(log_dir)) {
    cat('\n', log_dir, ' -> ', 'No the directory, create it.')
    dir.create(log_dir)
  }

  log_df <- NULL
  for(i in seq_len(experiment_code_count)){
    experiment_code_i <- experiment_code[i]
    args <- c(experiment_code_i, mascot_xml_dir, mascot_txt_dir)  # Set args to vector
    allArgs <- c(path2script, args)  # Add python script path to parameters vector
    log_out <- tryCatch(
      {
        output <- system2(command, args = allArgs, stdout = TRUE) # R call python script by pass parameters vector
        cat('\n', i, '->', experiment_code_i, '->', 'success', '\n')
        c(experiment_code_i, 'success')
      },

      warning = function(w){ # process warning
        cat('\n', i, '->', experiment_code_i, '->', 'warning', '\n')
        print(w)
        log_i <- c(experiment_code_i, 'warning')
        return(log_i)
      },

      error = function(e){ # process error
        cat('\n', i, '->', experiment_code_i, '->', 'error', '\n')
        print(e)
        log_i <- c(experiment_code_i, 'error')
        return(log_i)
      }
    )
    log_df <- rbind(log_df, log_out)
  }

  colnames(log_df) <- c('Exp_no', 'Status')
  now_time <- Sys.time()
  now_time <- gsub(':', '-', now_time)
  log_df_file_name <- paste(now_time, 'log_of_extract_psites_score.txt')
  log_df_file_path <- normalizePath(file.path(log_dir, log_df_file_name), mustWork = FALSE)
  utils::write.table(log_df, log_df_file_path, sep = '\t', row.names = FALSE, quote = FALSE)

  cat('\n Program finish, please see result log to check status.', '->', log_df_file_path)

}
