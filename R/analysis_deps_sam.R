#' Differential expression analysis using SAM
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param group A factor representing groups.
#' @param nperms Number of permutations used to estimate false discovery rates.
#' @param log2_label A boolean value for representing whether or not the value is logarithmic, the default is FALSE.
#' @param rand if specified, the random number generator will be put in a reproducible state.
#' @param minFDR A numeric value for filtering significant genes, the default is 0.05.
#' @param samr_plot A boolean value for representing whether or not samr graph is plotted.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references R. Tibshirani, G. Chu, T. Hastie and Balasubramanian Narasimhan (2010). samr: SAM: Significance Analysis of Microarrays.\
#' Rpackage version 1.28. https://CRAN.R-project.org/package=samr
#'
#' @return A list containing results from sam analysis.
#' @export
#'
#' @examples
#' \dontrun{
#' sam_results_list = analysis_deps_sam(
#'   expr_data_frame, group, log2_label = FALSE,
#'   nperms = 100, rand = NULL, minFDR = 0.05,samr_plot = T
#' )
#' }

analysis_deps_sam <- function(expr_data_frame, group, log2_label = FALSE,
                              nperms = 100, rand = NULL, minFDR = 0.05,
                              samr_plot = TRUE){
  requireNamespace('samr')
  requireNamespace('stats')
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


  # construct the samr data
  sam_data <- list(x = expr_Valule, y = as.numeric(as.factor(group)),
                   geneid = expr_ID, genenames = expr_ID, logged2=TRUE)

  group_nlevels <- nlevels(group)
  if(group_nlevels < 2){
    cat('\n', 'Groups are less than one.', '\n')
    stop('')
  }

  if(group_nlevels == 2){
    resp_type <- "Two class unpaired"
  }else{
    resp_type <- "Multiclass"
  }
  cat('\n', resp_type, '\n')
  samr_obj <- samr::samr(sam_data, resp.type = resp_type, nperms = nperms, random.seed = rand)

  # Compute the delta values
  delta_table <- samr::samr.compute.delta.table(samr_obj)

  # Determine a FDR cut-off
  index_less_than_min_FDR <- which(delta_table[,5] < minFDR)
  if(length(index_less_than_min_FDR) < 1){
    cat('\n', 'Not found appropiate cutoff less than specific minimum FDR.')
    stop('')
  }else{
    delta_index <- index_less_than_min_FDR[1]
    delta <- delta_table[delta_index,1]
  }


  if(samr_plot){
    cat('\n', 'Plot samr plot to view DEPs (or DEGs) distribution.')
    samr::samr.plot(samr_obj, delta)
  }

  # Extract significant genes at the cut-off delta
  siggenes_table <- samr::samr.compute.siggenes.table(samr_obj, delta, sam_data, delta_table, all.genes = F)
  genes_up_n <- siggenes_table$ngenes.up
  if(genes_up_n > 0){
    genes_up_df <- data.frame(siggenes_table$genes.up)
    genes_up_df_col <- ncol(genes_up_df)
    genes_up_df <- genes_up_df[,c(3,7:genes_up_df_col)]
    genes_up_df_col <- ncol(genes_up_df)
    genes_up_df[,genes_up_df_col] <- as.numeric(genes_up_df[,genes_up_df_col])/100
    genes_up_df_colnames <- colnames(genes_up_df)
    colnames(genes_up_df) <- c('ID', genes_up_df_colnames[-c(1,genes_up_df_col)], 'qvalue')

  }else{
    genes_up_df <- NULL
  }

  genes_lo_n <- siggenes_table$ngenes.lo
  if(genes_lo_n > 0){
    genes_lo_df <- data.frame(siggenes_table$genes.lo)
    genes_lo_df_col <- ncol(genes_lo_df)
    genes_lo_df <- genes_lo_df[,c(3,7:genes_lo_df_col)]
    genes_lo_df_col <- ncol(genes_lo_df)
    genes_lo_df[,genes_lo_df_col] <- as.numeric(genes_lo_df[,genes_lo_df_col])/100
    genes_lo_df_colnames <- colnames(genes_lo_df)
    colnames(genes_lo_df) <- c('ID', genes_lo_df_colnames[-c(1,genes_lo_df_col)], 'qvalue')
  }else{
    genes_lo_df <- NULL
  }

  sam_result_list <- list(
    genes_up_df <- genes_up_df,
    genes_down_df <- genes_lo_df
  )

  return(sam_result_list)
}
