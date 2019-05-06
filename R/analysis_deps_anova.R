#' Differential expression analysis using ANOVA
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param group A factor representing groups.
#' @param log2_label A boolean value for representing whether or not the value is logarithmic, the default is FALSE.
#' @param return_padjust A boolean value for representing whether or not the pvalue is adjusted, the default is TRUE.
#' @param adjust_method Correction method, such as "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr".
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame containing ID, log2(FC) and pvalue.
#' @export
#'
#' @examples
#' \dontrun{
#' anova_result <- analysis_deps_anova(
#'   expr_data_frame, group, log2_label = FALSE,
#'   return_padjust = TRUE, adjust_method = 'BH'
#' )
#' }



analysis_deps_anova <- function(expr_data_frame, group, log2_label = FALSE, return_padjust = TRUE, adjust_method = 'BH'){
  requireNamespace('stats')

  group_nlevels <- nlevels(group)
  if(group_nlevels < 2){
    cat('\n', 'Not found pairwise comparison.', '\n')
    stop('')
  }

  expr_ID <- as.vector(expr_data_frame[,1])
  if(!log2_label){
    expr_Valule <- log2(expr_data_frame[,-1]) # have to log
  }
  expr_Valule_row_duplicated <- apply(expr_Valule, 1, function(x){
    stats::var(x)
  })
  expr_Valule_col <- ncol(expr_Valule)
  duplicated_row_index <- which(expr_Valule_row_duplicated == 0)
  if(length(duplicated_row_index)>0){
    expr_ID <- expr_ID[-duplicated_row_index]
    expr_Valule <- expr_Valule[-duplicated_row_index,]
  }

  pvalue <- apply(expr_Valule, 1, function(y, x){
    y <- as.vector(unlist(y))
    fit <- stats::aov(y ~ x)
    fit.summary <- summary(fit)
    p <- fit.summary[[1]][1,5]
    p
  }, x = group)
  if(return_padjust){
    pvalue <- stats::p.adjust(pvalue, method <- adjust_method)
  }

  logFC <- apply(expr_Valule, 1, function(y, x){
      y <- as.vector(unlist(y))
      y_x_m <- tapply(y, x, mean)
      fc <- max(y_x_m) - min(y_x_m)
      fc
  }, x = group)


  anova_df <- data.frame(expr_ID, logFC, pvalue)
  colnames(anova_df) <- c('ID', 'logFC', 'pvalue')

  return(anova_df)
}
