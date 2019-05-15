#' Plot sequence logo based on list that consist of motifs and sequences.
#'
#' @param base_dir A path used for saving plots.
#' @param foreground_sequences_mapped_to_motifs A list that consist of motifs and sequences.
#' @param plot_min_seqs A numeric value for cutoff, sequences of motifs greater than the cutoff are plotted, the default is 5.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references (1) Omar Wagih (2017). ggseqlogo: A 'ggplot2' Extension for Drawing Publication-Ready Sequence Logos. R package version 0.1.\
#' https://github.com/omarwagih/ggseqlogo; (2) Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. \
#' R package version 1.3.0. https://CRAN.R-project.org/package=stringr
#'
#' @return Plot sequence logo based on list that consist of motifs and sequences.
#' The results will be saved in a folder named PhosMap_ggseqlogo in the BASE_DIR
#' parameter specified directory.
#'
#' @export
#' @examples
#' ftp_url <- "ftp://111.198.139.72:4000/pub/PhosMap_datasets/function_demo_data/foreground_sequences_mapped_to_motifs.RData"
#' load_data <- load_data_with_ftp(ftp_url, 'RData')
#' writeBin(load_data, "foreground_sequences_mapped_to_motifs.RData")
#' load("foreground_sequences_mapped_to_motifs.RData")
#'
#' BASE_DIR = getwd() # current working directory
#' BASE_DIR = normalizePath(BASE_DIR)
#' plot_seqlogo(BASE_DIR, foreground_sequences_mapped_to_motifs, plot_min_seqs = 25)
#'

plot_seqlogo <- function(base_dir, foreground_sequences_mapped_to_motifs, plot_min_seqs = 5){
  requireNamespace('ggseqlogo')
  requireNamespace('stringr')
  requireNamespace('grDevices')
  motifs_names <- names(foreground_sequences_mapped_to_motifs)
  motifs_names_count <- length(motifs_names)
  if(motifs_names_count < 1){
    cat('No motif results', '\n')
    stop('')
  }else{
    cat('Total motifs:', motifs_names_count, '\n')
  }

  ggseqlogo_dir <- normalizePath(file.path(base_dir, 'PhosMap_ggseqlogo'), mustWork = FALSE)
  if(!file.exists(ggseqlogo_dir)){
    dir.create(ggseqlogo_dir)
  }
  current_time <- Sys.time()
  current_time <- stringr::str_replace_all(current_time, ':', '-')
  current_time_dir <- normalizePath(file.path(ggseqlogo_dir, current_time), mustWork = FALSE)
  cat('Create a foder:', current_time_dir, '\n')
  dir.create(current_time_dir)

  plot_index <- 0
  for(i in seq_len(motifs_names_count)){
    motif_name <- motifs_names[i]
    motif_seq <- foreground_sequences_mapped_to_motifs[[i]]
    motif_seq_count <- length(motif_seq)
    if(motif_seq_count >= plot_min_seqs){
      motif_file_name <- paste(motif_seq_count, motif_name, ' seqs = ', motif_seq_count, '.pdf')
      motif_file_path <- normalizePath(file.path(ggseqlogo_dir, current_time, motif_file_name), mustWork = FALSE)

      grDevices::pdf(motif_file_path, height = 6, width = 6)
      print(ggseqlogo::ggseqlogo(motif_seq))
      grDevices::dev.off()

      plot_index <- plot_index + 1
      cat('Plot ', plot_index, ':', motif_name, '.\n', sep = '')

    }else{
      cat('Not Plot', motif_name, 'because the sequences of motif are less than', plot_min_seqs, '.\n')
    }
  }

  cat('\n', 'All sequence logoes are saved to', '\n  --', current_time_dir, '\n')
}
