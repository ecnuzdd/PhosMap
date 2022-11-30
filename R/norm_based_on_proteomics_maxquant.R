#' Normalizaiton for phosphoproteomics data from MaxQuant based on proteomics data
#' 
#' @param summary_phos_norm A data frame containing information required for all analysis.
#' @param proteomics_data A data frame of proteinGroups.txt.
#' @param experiment_code_file_path Experiment code file path for phosphoproteomics.
#' @param proteomics_experiment_file_path Experiment code file path for proteomics.
#' @param intensity_type Intensity type. The default is 'Intensity' and the options are 'iBAQ' and 'LFQ.intensity',depending on proteinGroups.txt.
#' @param min_unique_peptide Threshold for MaxQuant unique peptide[proteinGroups.txt]. The default is 1.
#' @param max_na_num Threshold for the number of missing values[proteinGroups.txt]. The default is 2.
#' @param norm_method Normalizaiton method[proteinGroups.txt]. The default is 'global' and the options are 'global' and 'median'.
#' @param impute_method Imputation method[proteinGroups.txt]. The default is 'minimum/10',the options are '0', 'minimum' and 'minimum/10'.
#' @return A result list. Elements are a data frame containing information required for all analysis and a pre-processed proteomics data.
#' @import utils stats
#' @export
#' @examples
#' 
#' ## Read phosphoproteomics data
#' \dontrun{
#' rawdata <- read.csv("Phospho (STY)Sites.txt",header=T,sep='\t')
#' 
#' ## Quality control for phosphoproteomics data 
#' qc_results <- qc_maxquant(rawdata, "./experiment_code_file.txt", min_score = 40, min_loc_prob = 0.75, max_na_num = 2)
#' qc_result <- qc_results[[1]]
#' qc_result_for_motifanalysis <- qc_results[[2]]
#' 
#' ## Normalizaiton, imputation and filtering
#' summary_phos_norm <- norm_maxquant(qc_result, qc_result_for_motifanalysis, norm_method = "global", impute_method = "minimum/10", top = 0.9)
#' 
#' ## Read proteomics data
#' proteomics_data <- read.csv("./proteinGroups.txt", sep = "\t")
#' 
#' results <- norm_based_on_proteomics(summary_phos_norm, proteomics_data, "./phosphorylation_exp_design_info.txt", "./phosphorylation_exp_design_info.txt")
#' summary_phos_norm_based_on_pro <- results[[1]]
#' pro_norm <- results[[2]]
#' }

norm_based_on_proteomics_maxquant <- function(
    summary_phos_norm,
    proteomics_data,
    experiment_code_file_path,
    proteomics_experiment_file_path,
    intensity_type = "Intensity",
    min_unique_peptide = 1,
    max_na_num = 2,
    norm_method = "global",
    impute_method = "minimum/10"
) {
  rawdata <- proteomics_data
  filenames <- read.csv(proteomics_experiment_file_path,header = T,sep='\t')
  
  rawdata <- rawdata[-which(rawdata$Reverse=="+"),]
  rawdata <- rawdata[-which(rawdata$Potential.contaminant=="+"),]
  rawdata <- rawdata[-which(rawdata$Protein.names == "" | rawdata$Gene.names == ""),]
  
  Genenames <- apply(data.frame(rawdata$Gene.names), 1, function(x){
    x <- strsplit(x, split = ';')[[1]]
    x[1]
  })
  
  varnames <- c()
  for(i in 1:nrow(filenames)){
    varnames[i] <- paste(intensity_type,filenames$Experiment_Code[i],sep = '.')
  }
  newdata <- rawdata[c('Unique.peptides', varnames)]
  newdata <- data.frame(Genenames, newdata)
  newdata2 <- newdata[which(newdata$Unique.peptides > min_unique_peptide),]
  NAnumthresig <- c()
  
  for (raw in 1:nrow(newdata2)) {
    NAnumthresig[raw] <- (sum(newdata2[raw,][-c(1,2)] == 0) <= max_na_num)
  }
  newdata3 <- newdata2[NAnumthresig,]
  
  if(norm_method == "global") {
    newdata4 <- sweep(newdata3[c(-1,-2)],2,apply(newdata3[c(-1,-2)],2,sum,na.rm=T),FUN="/")
  }
  if(norm_method == "median") {
    newdata4 <- sweep(newdata3[c(-1,-2)],2,apply(newdata3[c(-1,-2)],2,median,na.rm=T),FUN="/")
  }
  newdata4 <- newdata4 * 1e5
  
  newdata4[newdata4==0]<-NA
  df <- df1 <- newdata4
  if(impute_method=="0"){
    df[is.na(df)]<-0
  }else if(impute_method=="minimum"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)
  }else if(impute_method=="minimum/10"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)/10
  }
  df2 <- data.frame(newdata3$Genenames, df)
  colnames(df2) <- c("Symbol", filenames$Experiment_Code)
  
  phospho_data_topX <- summary_phos_norm[, c(-2, -3, -4)]
  colnames(phospho_data_topX)[1] <- "ID"
  
  phospho_data_topX_for_motifanalysis <- summary_phos_norm[, -1]
  
  data_frame_normalization_with_control_no_pair = normalize_phos_data_to_profiling(phospho_data_topX, df2,
                                                                                   experiment_code_file_path, proteomics_experiment_file_path,
                                                                                   pair_flag = FALSE)
  
  data_frame_normalization_with_control_no_pair_for_motifanalysis <- cbind(phospho_data_topX_for_motifanalysis[c(1,2,3)], data_frame_normalization_with_control_no_pair[-1])
  
  summarydf <- data.frame(data_frame_normalization_with_control_no_pair$ID, data_frame_normalization_with_control_no_pair_for_motifanalysis)
  colnames(summarydf) <- c("Position", colnames(data_frame_normalization_with_control_no_pair_for_motifanalysis))
  rownames(summarydf) <- rownames(data_frame_normalization_with_control_no_pair)
  results <- list(summarydf, df2)
  return(results)
}
