#' A simple PCA plot.
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param main The main title of plot.
#' @param point_cex a numerical value for point size.
#' @param point_col a color code or name for point color.
#' @param point_type point type, see points.
#' @param text_cex a numerical value for text size.
#' @export
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#'
#' @examples
#' \dontrun{
#' visualization_with_simple_pca(expr_data_frame, main = 'Simple PCA',
#' point_cex = 2, point_col = 'firebrick', point_type = 20, text_cex = 1)
#' }

visualization_with_simple_pca <- function(expr_data_frame,
                                          main = 'Simple PCA',
                                          point_cex = 2, point_col = 'firebrick', point_type = 20,
                                          text_cex = 1){
  requireNamespace('stats')
  requireNamespace('graphics')
  expr_ID <- as.vector(expr_data_frame[,1])
  expr_Valule <- log2(expr_data_frame[,-1]) # have to log
  testDat <- t(expr_Valule) # row -> sample, col -> variable
  pca <- stats::prcomp(((testDat)), center = TRUE, scale = TRUE)
  stats::screeplot(pca, type="lines")
  importance <- summary(pca)$importance
  PC1 <- importance[2,1]
  PC2 <- importance[2,2]
  PC1 <- round(PC1, 4)*100
  PC2 <- round(PC2, 4)*100

  pca_predict <- stats::predict(pca)
  pca_predict_2d <- pca_predict[,c(1,2)]
  ExpNames <- colnames(expr_Valule)
  rownames(pca_predict_2d) <- ExpNames

  #background
  xlim <- c(floor(min(pca_predict_2d[,1]))-5, ceiling(max(pca_predict_2d[,1]))+5)
  ylim <- c(floor(min(pca_predict_2d[,2]))-5, ceiling(max(pca_predict_2d[,2]))+5)
  xlab <- paste("PC1 (", PC1, "%)", sep = "")
  ylab <- paste("PC2 (", PC2, "%)", sep = "")
  graphics::plot(pca_predict_2d, type = "n", xlim = xlim, ylim = ylim, lwd = 2, xlab = xlab, ylab = ylab, main = main)
  graphics::points(pca_predict_2d, pch = point_type, col = point_col, cex = point_cex)
  graphics::text(pca_predict_2d, ExpNames, pos = 4, cex = text_cex)
}
