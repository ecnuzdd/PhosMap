#' Convert data frame of motif to the sequence pattern
#'
#' @param foreground_sequence A vector for AA sequences with fixed length as foreground input.
#' @param background_sequence A vector for AA sequences with fixed length as background input.
#' @param min_sequence_count A numeric for the minimum sequence number assigned to a motif.
#' @param min_pvalue A numeric for the minimum pvalue for found motif.
#' @param center A character for center of k-mer.
#' @param width A numeric for specific k-mer.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @references Omar Wagih (2014). rmotifx: An iterative statistical approach to the discovery of biological sequence motifs. R package version 1.0.
#'
#'
#' @return A list for information summary of searching mortif
#' @export
#'
#' @examples
#' \dontrun{
#' seach_motif_pattern(
#'     foreground_sequence,,
#'     background_sequence,
#'     min_sequence_count = 1,
#'     min_pvalue = 0.01,
#'     center = 'S',
#'     width = 15
#' )
#' }

seach_motif_pattern <- function(
  foreground_sequence,
  background_sequence,
  min_sequence_count = 1,
  min_pvalue = 0.01,
  center = 'S',
  width
){

  requireNamespace('stats')
  AA_LIST = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

  raw_foreground_sequence <- foreground_sequence
  raw_background_sequence <- background_sequence

  foreground_size <- length(foreground_sequence)
  background_size <- length(background_sequence)

  center_index <- ceiling(width/2)

  motif_coordinate <- matrix(0, 0, 2)
  pvalue <- c()

  # if there are no sequences in updated foreground_sequence or background_sequence, exit loop
  while(length(foreground_sequence) != 0 & length(background_sequence) != 0){
    # Construct foreground_pwm and background_pwm
    foreground_pwm <- construct_pwm(
      foreground_sequence,
      width,
      frequency_flag = TRUE
    )
    background_pwm <- construct_pwm(
      background_sequence,
      width,
      frequency_flag = FALSE # using probability distribution for binomial function
    )

    # Compute occurrence probability distribution of aa pattern in foreground using binomial function
    binomial_matrix <- 1 - stats::pbinom(foreground_pwm-1, length(foreground_sequence), background_pwm, lower.tail=TRUE)

    # replace 0 probability to 1e-16
    binomial_matrix[binomial_matrix == 0] <- 1e-16

    # flag central and determinated AAs and these AAs not meeting min_sequence_count
    binomial_matrix[,center_index] <- 1
    binomial_matrix[motif_coordinate] <- 1
    binomial_matrix[foreground_pwm < min_sequence_count] <- 1

    # search AAs equal to min_binomial_probability and less than min_pvalue
    min_binomial_probability <- min(binomial_matrix)
    min_binomial_probability_coordinate_matrix <- which(binomial_matrix == min_binomial_probability & binomial_matrix < min_pvalue, arr.ind=TRUE)

    # if there are no AAs meeting above reqirements, exit loop
    if(nrow(min_binomial_probability_coordinate_matrix) == 0){
      break
    }

    # search best match
    if(nrow(min_binomial_probability_coordinate_matrix) >= 1){
      max_weight_index <- which.max(foreground_pwm[min_binomial_probability_coordinate_matrix])
      aa_flags <- rownames(min_binomial_probability_coordinate_matrix)[max_weight_index]
      # update min_binomial_probability_coordinate_matrix
      min_binomial_probability_coordinate_matrix <- min_binomial_probability_coordinate_matrix[max_weight_index,]
      min_binomial_probability_coordinate_matrix <- t(as.matrix(min_binomial_probability_coordinate_matrix))
      rownames(min_binomial_probability_coordinate_matrix) <- aa_flags
    }
    row_index <- min_binomial_probability_coordinate_matrix[1,1]
    col_index <- min_binomial_probability_coordinate_matrix[1,2]
    # update foreground_sequence and background_sequence by extracting meeting min_pvalue
    foreground_sequence <- foreground_sequence[which(substr(foreground_sequence, col_index, col_index) == aa_flags)]
    background_sequence <- background_sequence[which(substr(background_sequence, col_index, col_index) == aa_flags)]

    # store specific aa index for motif pattern
    motif_coordinate <- rbind(motif_coordinate, c(row_index, col_index))
    # cat('\n', length(foreground_sequence), length(background_sequence), '\n')

    pvalue <- c(pvalue, min_binomial_probability)
  }

  # Motif data: data frame with amino acid and positions for the motif
  motif_coordinate_data_frame <- as.data.frame(motif_coordinate)
  colnames(motif_coordinate_data_frame) <- c('aa', 'index')


  motif_coordinate_data_frame$aa <- AA_LIST[motif_coordinate_data_frame$aa]

  # Not found motif
  if(nrow(motif_coordinate_data_frame) == 0){
    return(NULL)
  }

  # Convert data frame of motif to the sequence pattern
  motif_pattern <- motif_data_frame_to_sequence(motif_coordinate_data_frame, center, width)

  # get the motif score
  if(length(pvalue)>0){
    motif_pattern_score <-  sum(-log10(pvalue))
  }else{
    motif_pattern_score <- NA
  }

  # count matches in foreground and background, respectively.
  foreground_matches <- length(grep(motif_pattern, raw_foreground_sequence))
  background_matches <- length(grep(motif_pattern, raw_background_sequence))


  return(
    list(
      foreground = foreground_sequence,
      motif_coordinate_data_frame = motif_coordinate_data_frame,
      motif_pattern_score = motif_pattern_score,
      motif_pattern = motif_pattern,
      foreground_matches = foreground_matches,
      foreground_size = foreground_size,
      background_matches = background_matches,
      background_size = background_size
    )
  )

}



