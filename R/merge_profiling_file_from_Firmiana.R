#' Merge profiling files downloaded from Firmiana.
#'
#' Filter data based on US (Unique and Ionscore > 20) peptide counts then merge profiling files.
#'
#'
#' @param firmiana_gene_dir a folder containing gene identification results as input.
#' @param US_cutoff a numerical value as a cutoff to filter data, the default is 1.
#' @param experiment_gene_file_path a file path for storing experiemnt design of proteomics data.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A merged data frame after filtering (US_cutoff) and replacing NAs to zeros.
#' @export
#'
#' @examples
#' \dontrun{
#' merged_df = merge_profiling_file_from_Firmiana(firmiana_gene_dir, US_cutoff = 1)
#' }
merge_profiling_file_from_Firmiana <- function(firmiana_gene_dir, US_cutoff = 1, experiment_gene_file_path){
  requireNamespace('utils')
  DATA_DIR <- normalizePath(firmiana_gene_dir, mustWork = FALSE)
  if(!file.exists(DATA_DIR)){
    cat(DATA_DIR, ' -> ', 'No the file')
    stop('')
  }
  data_list <- list()
  file_names <- list.files(path = DATA_DIR, pattern = '.txt')
  file_names_count <- length(file_names)
  if(length(file_names_count)<1){
    stop('The directory of ', DATA_DIR, ' has no files.')
  }

  exp_names <- apply(data.frame(file_names), 1, function(x){
    x <- strsplit(x, split = '_')[[1]][1]
    x
  })

  experiment_code <- utils::read.table(experiment_gene_file_path, header = TRUE, sep = '\t', stringsAsFactors = NA)
  experiment_code <- as.vector(unlist(experiment_code$Experiment_Code))

  index_of_match <- match(experiment_code, exp_names)
  na_index <- which(is.na(index_of_match))
  na_count <- length(na_index)
  if(na_count > 0){
    na_experiment_code <- experiment_code[na_index]
    cat(
      '\n',
      na_experiment_code,
      'not in',
      DATA_DIR
    )
    stop('')
  }

  exp_names <- exp_names[index_of_match]
  file_names <- file_names[index_of_match]
  file_names_count <- length(file_names)


  # Table headers of input data
  # "Gene.ID" "Symbol" "Annotation"  "Modification" "Description"
  # "Protein.GI" "Protein.Num" "Area" "FoT.1e.6." "iBAQ"
  # "Peptide.Num" "Unique.Peptide.Num"  "Strict.Peptide.Num"  "US.Peptide.Num"  "Identified.Proteins.Num"
  # "Unique.Proteins.Num"

  # New table headers of input data
  file_data_colnames <- c(
    "Gene_ID", "Symbol", "Annotation", "Modification", "Description",
    "Protein_GI",  "Protein_Num", "Area", "FoT5", "iBAQ",
    "Peptide_Num", "UPeptide_Num",  "SPeptide_Num",  "USPeptide_Num",  "Identified_Proteins_Num", "Unique_Proteins_Num"
  )
  kept_colnames <- c(
    "Symbol", "iBAQ", "USPeptide_Num"
  )
  kept_colnames_index <- match(kept_colnames, file_data_colnames)
  cat('\n Merge profiling files downloaded from Firmiana.')
  cat('\n Total files: ', file_names_count)
  for(i in seq_len(file_names_count)){
    file_name <- file_names[i]
    file_path <- normalizePath(file.path(DATA_DIR, file_name))
    file_data <- utils::read.delim(file_path, header = TRUE, stringsAsFactors = NA, sep = '\t')
    colnames(file_data) <- file_data_colnames
    file_data <- file_data[, kept_colnames_index]

    index_of_US <- which(file_data$USPeptide_Num >= US_cutoff)
    file_data <- file_data[index_of_US, c(1,2)]
    exp_name <- exp_names[i]
    file_data_colnames.i <- colnames(file_data)
    file_data_colnames.i <- paste(exp_name, file_data_colnames.i, sep = '_')
    file_data_colnames.i[1] <- 'Symbol'
    colnames(file_data) <- file_data_colnames.i
    data_list[[i]] <- file_data
    cat('\n Read and filter: ', i, '/', file_names_count)
  }
  attr(data_list, 'names') <- exp_names

  data_list_count <- length(data_list)
  merge_df <- data_list[[1]]
  merge_df_colnames <- colnames(merge_df)
  cat('\n Merge_completed: ', 1, '/', data_list_count)
  if(data_list_count>1){
    for(i in 2:data_list_count){
      tmp_merge_df <- data_list[[i]]
      merge_df <- merge(merge_df, tmp_merge_df, by = 'Symbol', all = TRUE)
      cat('\n merge_complete: ', i, '/', data_list_count)
    }
  }
  Symbol <- as.vector(merge_df[,1])
  Value <- as.matrix(merge_df[,-1])
  index_of_NA <- which(is.na(Value))
  if(length(index_of_NA)>0){
    Value[index_of_NA] <- 0
  }
  colnames(Value) <- exp_names
  merge_df_no_NA <- data.frame(Symbol, Value)
  return(merge_df_no_NA)
}








