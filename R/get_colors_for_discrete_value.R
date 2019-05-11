#' Generate custom colors from discrete values for heatmaps.
#'
#'
#' @param color_intervals_list a list for building color intervals.
#' @param value_intervals_list a list for building value intervals.
#' @author Dongdong Zhan and Mengsha Tong
#'
#'
#' @return A vectors containing color distributions.
#' @export
#'
#' @examples
#' \dontrun{
#' value_intervals_list <- list(
#' seq(-4, -2, 0.2),
#' seq(-2, -1, 0.2),
#' seq(-1, 1, 0.2),
#' seq(1, 2, 0.2),
#' seq(2, 4, 0.2)
#' )
#' color_intervals_list <- list(
#'   c('blue', '#33CCFF'),
#'   c('#33CCFF', 'green'),
#'   c('green', 'white', '#FF6600'),
#'   c('#FF6600', 'red'),
#'   c('red', 'firebrick')
#' )
#' colors <- get_colors_for_discrete_value(
#'   color_intervals_list,
#'   value_intervals_list
#' )
#' }
#'

get_colors_for_discrete_value <- function(color_intervals_list, value_intervals_list){
  requireNamespace('grDevices')
  color_intervals_list_len <- length(color_intervals_list)
  value_intervals_list_len <- length(value_intervals_list)
  if(color_intervals_list_len != value_intervals_list_len){
    cat('\n', 'The intervals of colors and values are not equal.', '\n')
    stop('')
  }else{
    colors <- NULL
    for(i in seq_len(color_intervals_list_len)){
      breaks_i <- value_intervals_list[[i]]
      colors_i <- color_intervals_list[[i]]
      if(i == 1){
        colors_v <- grDevices::colorRampPalette(colors_i)(length(breaks_i)-1)
      }else{
        colors_v <- grDevices::colorRampPalette(colors_i)(length(breaks_i))

      }
      # cat('\n', colors_v, '\n')
      colors <- c(colors, colors_v)
    }
    colors <- colors[which(!duplicated(colors))]
    return(colors)
  }
}
