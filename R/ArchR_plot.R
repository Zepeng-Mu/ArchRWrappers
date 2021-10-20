###############################################################################
## Script purpose: Helper functions for plotting ArchR data
## Author: Zepeng Mu
###############################################################################

#' Plot Reduced Dimensions
#'
#' @param mtrx 
#' @param maxDim 
#' @param color 
#' @param discrete 
#'
#' @return
#' @export
#'
#' @examples
plotRedDim <- function(mtrx, maxDim = 30, color = NULL, discrete = F) {
  nCol <- ceiling(sqrt(maxDim))
  
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
