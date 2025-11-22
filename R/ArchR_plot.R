###############################################################################
## Script purpose: Helper functions for plotting ArchR data
## Author: Zepeng Mu
###############################################################################

#' A minimal ggplot2 theme for publication-quality plots
#'
#' This function returns a clean, classic ggplot2 theme with customized text sizes
#' suitable for scientific publications. It removes the legend title and sets
#' consistent font sizes for axis text and titles.
#'
#' @param ... Additional arguments passed to \code{theme_classic}
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   my_theme()
#' }
my_theme <- function(...) {
  theme_classic(base_size = 7) +
    theme(
      axis.text = element_text(size = unit(7, "pt")),
      axis.title = element_text(size = unit(7, "pt")),
      legend.title = element_blank(),
      legend.text = element_text(size = unit(7, "pt"))
    )
}

#' Plot reduced dimensions in a grid layout
#'
#' This function creates a grid of scatter plots showing consecutive pairs of
#' reduced dimensions (e.g., LSI1 vs LSI2, LSI2 vs LSI3, etc.), useful for
#' visualizing the structure of dimensionality reduction results.
#'
#' @param mtrx A numeric matrix containing reduced dimensions, with cells as rows
#'   and dimensions as columns
#' @param maxDim Integer specifying the maximum number of dimensions to plot.
#'   This cannot be larger than the number of columns in mtrx
#' @param color A vector with length equal to the number of rows in mtrx, used
#'   for coloring each point. Default is NULL
#' @param discrete Logical indicating whether the color variable should be treated
#'   as discrete. Default is TRUE
#' @param nCol Integer specifying the number of columns in the plot grid. If NULL,
#'   automatically calculated as ceiling(sqrt(maxDim))
#'
#' @return A patchwork object containing the grid of dimension plots
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot first 20 LSI dimensions
#' lsiMat <- getReducedDims(projHeme, reducedDims = "IterativeLSI")
#' plotRedDim(
#'   mtrx = lsiMat,
#'   maxDim = 20,
#'   color = projHeme$Clusters,
#'   discrete = TRUE
#' )
#' }
plotRedDim <- function(
  mtrx,
  maxDim = 30,
  color = NULL,
  discrete = T,
  nCol = NULL
) {
  if (maxDim > ncol(mtrx)) {
    stop(
      "maxDim cannot be larger be larger than the number of factors (column numbers) in mtrx"
    )
  }

  if (is.null(nCol)) {
    nCol <- ceiling(sqrt(maxDim))
  }

  plotList <- ArchR:::.safelapply(
    1:(maxDim - 1),
    function(x) {
      ggPoint(
        x = mtrx[, x],
        y = mtrx[, x + 1],
        xlabel = "Dim" %&% x,
        ylabel = "Dim" %&% (x + 1),
        size = 0.1,
        color = color,
        discreteSet = "kelly",
        baseSize = 7,
        extend = 0.01,
        discrete = discrete
      ) +
        theme(
          plot.margin = unit(c(0, 0.15, 0, 0), "in"),
          aspect.ratio = 1
        )
    }
  )

  outPlot <- Reduce("+", plotList) +
    plot_layout(ncol = nCol, guides = "collect") &
    theme(legend.position = "bottom")

  return(outPlot)
}
