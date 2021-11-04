###############################################################################
## Script purpose: Helper functions for plotting ArchR data
## Author: Zepeng Mu
###############################################################################

#' A ggplot2 theme
my_theme <- function(...) {
  theme_classic(base_size = 7) +
    theme(
      axis.text = element_text(size = unit(7, "pt")),
      axis.title = element_text(size = unit(7, "pt")),
      legend.title = element_blank(),
      legend.text = element_text(size = unit(7, "pt"))
    )
}

#' Plot Reduced Dimensions
#'
#' @param mtrx 
#' @param maxDim The maximum number of dimensions to plot. This cannot be larger than the number of factors in mtrx
#' @param color A vector with equal length as row number of mtrx for coloring each point
#' @param discrete Whether the color should be discrete. Boolean. Default is TRUE
#' @param nCol 
#'
#' @return
#' @export
#'
#' @examples
plotRedDim <- function(mtrx, maxDim = 30, color = NULL, discrete = T, nCol = NULL) {
  if (maxDim > ncol(mtrx)) {
    stop("maxDim cannot be larger be larger than the number of factors (column numbers) in mtrx")
  }
  
  if (is.null(nCol)) {
    nCol <- ceiling(sqrt(maxDim))
  }

  plotList <- ArchR:::.safelapply(
    1:(maxDim - 1),
    function(x) {
      ggPoint(x = mtrx[, x],
              y = mtrx[, x + 1],
              xlabel = "Dim"%&%x,
              ylabel = "Dim"%&%(x + 1),
              size = 0.1,
              color = color,
              discreteSet = "kelly",
              baseSize = 7,
              extend = 0.01,
              discrete = discrete) +
        theme(plot.margin = unit(c(0, 0.15, 0, 0), "in"),
              aspect.ratio = 1)
    })
  
  outPlot <- Reduce("+", plotList) +
    plot_layout(ncol = nCol, guides = "collect") &
    theme(legend.position = "bottom")
  
  return(outPlot)
}
