#' A simple t-SNE plot.
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param group A factor for group information.
#' @param main The main title of plot.
#' @param perplexity A numerical value for perplexity, the default is 10.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @export
#'
#'
#' @examples
#' \dontrun{
#' visualization_with_simple_tsne(
#'   expr_data_frame,
#'   group,
#'   main = 'Simple t-SNE',
#'   perplexity = 10
#' )
#' }

visualization_with_simple_tsne <- function(expr_data_frame, group, main = 'Simple t-SNE', perplexity = 10)
{
  requireNamespace('Rtsne')
  requireNamespace('grDevices')
  requireNamespace('graphics')

  expr_ID <- as.vector(expr_data_frame[,1])
  expr_Valule <- log2(expr_data_frame[,-1]) # have to log
  rtsne_df <- data.frame(group, t(expr_Valule))
  rtsne_df_unique <- unique(rtsne_df)
  colors <- grDevices::rainbow(length(unique(group)))
  names(colors) <- unique(group)
  # Executing the algorithm on curated data
  system.time({
    tsne <- Rtsne::Rtsne(rtsne_df_unique[,-1], dims = 2, perplexity = perplexity, verbose = TRUE, max_iter = 500)
  })
  graphics::plot(tsne$Y, t = 'n', main = main, xlab = 'Dim 1', ylab = 'Dim 2')
  graphics::text(tsne$Y, labels = group, col = colors[group])
}
