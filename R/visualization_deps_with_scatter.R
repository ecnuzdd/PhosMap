#' Visualize differentially expressed results with scatter
#'
#' @param deps_data a data frame containing ID, logFC and pvalue.
#' @param minFC a numeric for the minimum fold change.
#' @param minPvalue a numeric for the significance cutoff.
#' @param main an overall title for the plot.
#' @param show_text a boolean value representing whether or not the text is showed, the default is FALSE.
#' @param min_up_text cutoff value for showing up-IDs. Only IDs with lower than min_up_text are showed.
#' @param min_down_text cutoff value for showing down-IDs. Only IDs with lower than min_down_text are showed.
#'
#' @author Dongdong Zhan and Mengsha Tong
#' @export
#'
#'
#' @return A scatter plot for showing differentially expressed results.
#'
#' @examples
#' \dontrun{
#' visualization_deps_with_scatter(
#'   deps_data,
#'   minFC = 2,
#'   minPvalue = 0.05,
#'   main = 'Differentially expressed proteins',
#'   show_text = FALSE,
#'   min_up_text = 15,
#'   min_down_text = 15
#' )
#' }

visualization_deps_with_scatter <- function(
  deps_data,
  minFC = 2,
  minPvalue = 0.05,
  main = 'Differentially expressed proteins',
  show_text = FALSE,
  min_up_text = 15,
  min_down_text = 15
){
  requireNamespace('graphics')
  requireNamespace('stats')
  x_v <- deps_data$logFC
  x_v_max <- max(x_v)
  x_v_right <- ceiling(x_v_max)
  x_v_min <- min(x_v)
  x_v_left <- floor(x_v_min)

  x_up <- log2(minFC)
  x_down <- log2(1/minFC)

  zero_index <- which(deps_data$pvalue==0)
  zero_index_count <- length(zero_index)
  if(zero_index_count){
    minimum_p <- min(deps_data$pvalue[-zero_index])
    min <- minimum_p/10
    max <- minimum_p-minimum_p/10
    minimum_p_new <- stats::runif(zero_index_count, min = min, max = max)
    deps_data$pvalue[zero_index] <- minimum_p_new
  }

  y_v <- (-log10(deps_data$pvalue))
  y_v_max <- max(y_v)
  y_v_up <- ceiling(y_v_max)
  y_v_sig <- (-log10(minPvalue))


  index_of_up <- which(x_v > x_up & y_v > y_v_sig)
  index_of_down <- which(x_v < x_down & y_v > y_v_sig)


  graphics::plot(x_v, y_v,
       xlim = c(x_v_left, x_v_right), ylim = c(0, y_v_up),
       xlab = 'log2(FC)', ylab = '-log10(pvalue)', main = main)
  graphics::abline(h = y_v_sig, lty = 'dotdash', col = 'firebrick', lwd = 2)
  graphics::abline(v = x_up, lty = 'dotdash', col = 'firebrick', lwd = 2)
  graphics::abline(v = x_down, lty = 'dotdash', col = 'firebrick', lwd = 2)

  graphics::points(x_v[index_of_up], y_v[index_of_up], pch = 20, col = 'red')
  graphics::points(x_v[index_of_down], y_v[index_of_down], pch = 20, col = 'blue')

  if(show_text){
    s <- as.vector(deps_data$ID)
    s_up <- s[index_of_up]
    x_v_up_set <- x_v[index_of_up]
    x_v_up_set_order <- order(x_v_up_set, decreasing = TRUE)
    y_v_up_set <- y_v[index_of_up]
    y_v_up_set_order <- order(y_v_up_set, decreasing = TRUE)

    index_up_set <- intersect(x_v_up_set_order[seq_len(min_up_text)], y_v_up_set_order[seq_len(min_up_text)])
    graphics::text(x_v_up_set[index_up_set], y_v_up_set[index_up_set], s_up[index_up_set], pos = 3, cex = 0.6)

    s_down <- s[index_of_down]
    x_v_down_set <- x_v[index_of_down]
    x_v_down_set_order <- order(x_v_down_set, decreasing = FALSE)
    y_v_down_set <- y_v[index_of_down]
    y_v_down_set_order <- order(y_v_down_set, decreasing = TRUE)

    index_down_set <- intersect(x_v_down_set_order[seq_len(min_down_text)], y_v_down_set_order[seq_len(min_down_text)])
    graphics::text(x_v_down_set[index_down_set], y_v_down_set[index_down_set], s_down[index_down_set], pos = 3, cex = 0.6)
  }
}
