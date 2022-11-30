#' Normalizaiton, imputation and filtering for phosphoproteomics data from MaxQuant.
#' 
#' @param qc_result A data frame containing quality control result.
#' @param qc_result_for_motifanalysis A data frame containing information required for motif analysis.
#' @param norm_method Normalizaiton method. The default is 'global' and the options are 'global' and 'median'.
#' @param impute_method Imputation method. The default is 'minimum/10' and the options are '0', 'minimum' and 'minimum/10'.
#' @param percent_of_kept_sites A numeric value representing a cutoff used for filter psites. The default is 3/4.
#' @import stats
#' @export
#' @return A data frame containing information required for all analysis.
#' @examples
#' \dontrun{
#' rawdata <- read.csv("Phospho (STY)Sites.txt",header=T,sep='\t')
#' 
#' ## Quality control for phosphoproteomics data 
#' qc_results <- qc_maxquant(rawdata, "./experiment_code_file.txt", min_score = 40, min_loc_prob = 0.75, max_na_num = 2)
#' qc_result <- qc_results[[1]]
#' qc_result_for_motifanalysis <- qc_results[[2]]
#' 
#' summary_phos_norm <- norm_maxquant(qc_result, qc_result_for_motifanalysis, norm_method = "global", impute_method = "minimum/10", top = 0.9)
#' }

norm_maxquant <- function(
    qc_result,
    qc_result_for_motifanalysis,
    norm_method = "global",
    impute_method = "minimum/10",
    percent_of_kept_sites = 3/4
    
) {
  newdata3out <- qc_result
  newdata3motif <- qc_result_for_motifanalysis
  top <- percent_of_kept_sites
  if(norm_method == "global") {
    newdata4 <- sweep(newdata3out[-1],2,apply(newdata3out[-1],2,sum,na.rm=T),FUN="/")
  }
  if(norm_method == "median") {
    newdata4 <- sweep(newdata3out[-1],2,apply(newdata3out[-1],2,median,na.rm=T),FUN="/")
  }
  newdata4 <- newdata4 *1e5
  newdata4[newdata4==0]<-NA
  df <- df1 <- newdata4
  if(impute_method=="0"){
    df[is.na(df)]<-0
  }else if(impute_method=="minimum"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)
  }else if(impute_method=="minimum/10"){
    df[is.na(df)]<-min(df1,na.rm = TRUE)/10
  }
  
  dfmotif <- data.frame(newdata3motif$AA_in_protein, newdata3motif$Sequence, newdata3motif$ID, df)
  colnames(dfmotif) <- colnames(newdata3motif)
  rownames(dfmotif) <- seq(nrow(dfmotif))
  df <- data.frame(newdata3out$ID, df)
  colnames(df) <- colnames(newdata3out)
  
  phospho_data_topX = keep_psites_with_max_in_topX(df, percent_of_kept_sites = top)
  percent_of_kept_sites_str <- paste('top', top*100, '%', sep = '')
  Value <- dfmotif[,-c(1,2,3)]
  Value_rowmax <- apply(Value, 1, function(x){
    x <- as.vector(unlist(x))
    max(x)
  })
  index_of_Value_rowmax_desc <- order(Value_rowmax, decreasing = TRUE)
  count_of_kept_sites <- round(nrow(Value)*top)
  index_of_Value_rowmax_desc_kept <- index_of_Value_rowmax_desc[seq_len(count_of_kept_sites)]
  phospho_data_topX_for_motifanalysis <- dfmotif[index_of_Value_rowmax_desc_kept,]
  summarydf <- data.frame(phospho_data_topX$ID, phospho_data_topX_for_motifanalysis)
  colnames(summarydf) <- c("Position", colnames(phospho_data_topX_for_motifanalysis))
  rownames(summarydf) <- rownames(phospho_data_topX)
  
  return(summarydf)
}
