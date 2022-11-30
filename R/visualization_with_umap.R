#' A umap plot.
#'
#' @param expr_data_frame A data frame containing ID and quantification value.
#' @param group A factor for group information.
#' @param main The main title of plot.
#' @param n_neighbors A numerical value for the size of local neighborhood, the default is 10.
#' @import uwot graphics
#' @export
#' 
#' @return A umap plot.

visualization_with_umap <- function(
    expr_data_frame,
    group,
    main = "UMAP",
    n_neighbors = 10
) {
  requireNamespace('uwot')
  requireNamespace('grDevices')
  
  expr_data_frame_2=t(expr_data_frame[,-1])
  expr_data_frame_2 <- as.data.frame(expr_data_frame_2)
  
  expr_data_group <- group
  expr_data_umap <- umap(expr_data_frame_2,n_neighbors = n_neighbors)
  color_group_unique <- grDevices::rainbow(length(unique(expr_data_group)))
  plot(expr_data_umap,col=rep(color_group_unique,table(expr_data_group)),pch=16,asp=1,xlab = "UMAP_1",ylab = "UMAP_2",main = main)
  abline(h=0,v=0,lty=2,col="gray")
  legend("topright",title = "Group",inset = 0.01,
         legend = unique(expr_data_group), pch=16,
         col = color_group_unique,
         bg='white'
  )
}