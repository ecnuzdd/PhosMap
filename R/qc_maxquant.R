#' Quality control for phosphoproteomics data from MaxQuant.
#' 
#' @param data_frame A data frame of Phospho (STY)Sites.txt.
#' @param experiment_code_file_path Experiment code file path.
#' @param min_score Threshold for MaxQuant score. The default is 40.
#' @param min_loc_prob Threshold for MaxQuant Localization.prob. The default is 0.75.
#' @param max_na_num Threshold for the number of missing values. The default is 2.
#' @import utils stringr
#' @export
#' @return A result list. Elements are a data frame containing quality control result and a data frame containing information required for motif analysis.
#' @examples
#' \dontrun{
#' rawdata <- read.csv("Phospho (STY)Sites.txt",header=T,sep='\t')
#' qc_results <- qc_maxquant(rawdata, "./experiment_code_file.txt", min_score = 40, min_loc_prob = 0.75, max_na_num = 2)
#' qc_result <- qc_results[[1]]
#' qc_result_for_motifanalysis <- qc_results[[2]]
#' }

qc_maxquant <- function(
    data_frame,
    experiment_code_file_path,
    min_score = 40,
    min_loc_prob = 0.75,
    max_na_num = 2
) {
  requireNamespace('utils')
  rawdata <- data_frame
  filenames <- read.csv(experiment_code_file_path,header = T,sep='\t')
  rawdata <- rawdata[-which(rawdata$Reverse=="+"),]
  rawdata <- rawdata[-which(rawdata$Potential.contaminant=="+"),]
  rawdata <- rawdata[-which(rawdata$Protein.names == "" | rawdata$Gene.names == ""),]
  
  varnames <- c()
  for(i in 1:nrow(filenames)){
    varnames[i] <- paste('Intensity',filenames$Experiment_Code[i],sep = '.')
  }
  varnames <- gsub('-','.',varnames)
  myvar <- c('Localization.prob','Score', 'Leading.proteins','Gene.names','Amino.acid', 'Positions', 'Phospho..STY..Probabilities', 'Position.in.peptide', varnames)
  newdata <- rawdata[myvar]
  
  Genenames <- apply(data.frame(newdata$Gene.names), 1, function(x){
    x <- strsplit(x, split = ';')[[1]]
    x[1]
  })
  Positions <- apply(data.frame(newdata$Positions), 1, function(x){
    x <- strsplit(x, split = ';')[[1]]
    x[1]
  })
  Leadingproteins <- apply(data.frame(newdata$Leading.proteins), 1, function(x){
    x <- strsplit(x, split = ';')[[1]]
    x[1]
  })
  newdata$Gene.names <- Genenames
  newdata$Positions <- Positions
  newdata$Leading.proteins <- Leadingproteins
  
  newdata2 <- newdata[which(newdata$Localization.prob >= min_loc_prob & newdata$Score >= min_score),]
  
  NAnumthresig <- c()
  for (row in 1:nrow(newdata2)) {
    NAnumthresig[row] <- (sum(newdata2[row,][-c(seq(8))] == 0) <= max_na_num)
  }
  newdata3 <- newdata2[NAnumthresig,]
  
  aminoposition <- paste0(newdata3$Amino.acid, newdata3$Positions)
  ID <- paste(newdata3$Gene.names, aminoposition, sep = '_')
  rowname <- paste(newdata3$Leading.proteins, ID, sep = '_')
  
  sequence <- gsub("\\(.*?\\)","",newdata3$Phospho..STY..Probabilities)
  Sequence <- c()
  for (i in seq(length(sequence))) {
    tmp <- unlist(str_split(sequence[i], pattern = ""))
    index <- newdata3$Position.in.peptide[i]
    tmp[index] <- tolower(tmp[index])
    
    Sequence[i] <- paste(tmp, collapse = "")
  }
  
  newdata3out <- newdata3[-c(seq(8))]  
  rownames(newdata3out) <- rowname
  
  newdata3motif <- data.frame(aminoposition, Sequence, newdata3$Leading.proteins, newdata3out)
  colnames(newdata3motif) <- c("AA_in_protein", "Sequence", "ID", filenames$Experiment_Code)
  newdata3out <- data.frame(ID, newdata3out)
  colnames(newdata3out) <- c("ID", filenames$Experiment_Code)
  newdata3motif <- newdata3motif[!duplicated(newdata3out$ID),]
  newdata3out <- newdata3out[!duplicated(newdata3out$ID),]
  result <- list(newdata3out, newdata3motif)
  return(result)
}
