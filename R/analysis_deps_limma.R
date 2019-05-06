#' Differential expression analysis using limma.
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param group A factor for representing groups.
#' @param comparison_factor A vector for comparison factor.
#' @param log2_label A boolean value for representing whether or not the value is logarithmic, the default is FALSE.
#' @param adjust_method method used to adjust the p-values for multiple testing. See p.adjust for the complete list of options, the default is "BH"
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression \
#' analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#'
#' @return A list containing results from limma analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' limma_results_list = analysis_deps_limma(
#'   expr_data_frame, group, comparison_statement,
#'   log2_label = FALSE, adjust_method = 'BH'
#' )
#' }
analysis_deps_limma <- function(expr_data_frame, group, comparison_factor,
                                log2_label = FALSE, adjust_method = 'BH'){
  requireNamespace('limma')
  requireNamespace('stats')
  # experiment_design_file_path = "D:\\Phosphate-data\\Bioinfomatics\\demo_data_from_WYN\\experiment_design_noPair.txt"
  # experiment_design_file = read.table(experiment_design_file_path, sep = '\t', header = T)
  # group = experiment_design_file$Group[experiment_design_file$Data_Type == 'Phospho']
  # group = paste('t', group, sep = '')
  # group = factor(group, levels = c('t0', 't10', 't30', 't120'))
  # expr_data_frame = data_frame_normalization_0

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
    #  Zero sample variances detected, have been offset away from zero
    expr_ID <- expr_ID[-duplicated_row_index]
    expr_Valule <- expr_Valule[-duplicated_row_index,]
  }
  # rownames(expr_Valule) <- expr_ID

  design <- stats::model.matrix(~ 0 + group)
  cat('\n', 'The matrix of experiment design.')
  print(design)
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(expr_Valule)
  # comparison_statement <- c('t10-t0', 't30-t0', 't120-t0')
  # comparison_statement <- c('t10-t0')
  group_levels <- comparison_factor
  group_levels_count <- length(group_levels)
  if(group_levels_count<2){
    cat('\n', 'Do not construct pairwise comparison pattern.')
    stop('')
  }else{
    comparison_statement <- NULL
    i_end <- group_levels_count - 1
    for(i in seq_len(i_end)){
      ctrl <- group_levels[i]
      j_start <- i + 1
      for(j in j_start:group_levels_count){
        treat <- group_levels[j]
        cs <- paste(treat, '-', ctrl, sep = '')
        comparison_statement <- c(comparison_statement, cs)
      }
    }
    cat('\n', 'The combination of pairwise comparison(s).')
    cat('\n', comparison_statement, '\n')
  }



  contrast.matrix <- limma::makeContrasts(contrasts = comparison_statement, levels = design)
  cat('\n', 'The matrix of comparison statement, compare other groups with control.')
  print(contrast.matrix) # the matrix of comparison statement, compare other groups with control.


  # step1
  fit <- limma::lmFit(expr_Valule, design)

  # step2
  fit2 <- limma::contrasts.fit(fit, contrast.matrix) # An important step.
  fit2 <- limma::eBayes(fit2)  # default no trend!


  # return(fit2)
  # step3
  alls <- limma::topTable(fit2, coef = 1, adjust.method = adjust_method, p.value = 1, number = Inf) # logFC = log(a/b) = log(a) - log(b) = A - B
  # results <- decideTests(fit2, method = "global", adjust.method = adjust_method, p.value = minPvalue, lfc = minFC)
  # vennDiagram(results)
  alls <- stats::na.omit(alls)

  # plot
  ID <- rownames(alls)
  logFC <- alls$logFC # log2
  pvalue <- alls$adj.P.Val

  result_df <- data.frame(ID, logFC, pvalue)

  return(result_df)
}




































