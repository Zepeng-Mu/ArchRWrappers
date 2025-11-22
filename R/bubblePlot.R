#' Create a bubble plot for visualizing marker features
#'
#' This function creates a bubble plot (also known as dot plot) to visualize marker
#' features across different groups. The size of each bubble represents the proportion
#' of cells expressing the feature, while the color represents the mean expression level.
#' Significant markers can be highlighted with black borders.
#'
#' @param propMtrx A numeric matrix containing the proportion of cells expressing each
#'   feature (rows) in each group (columns). Values should be between 0 and 1
#' @param meanMtrx A numeric matrix containing the mean expression values for each
#'   feature (rows) in each group (columns)
#' @param markerSet A SummarizedExperiment object containing marker statistics with
#'   "FDR" and "Log2FC" assays, typically from ArchR's getMarkerFeatures
#' @param useFeatures Character vector of feature names to include in the plot
#' @param labelFeatures Character vector of feature names to label in the plot.
#'   Default is the same as useFeatures
#' @param Log2FC Numeric threshold for log2 fold change to define significant markers.
#'   Default is 1
#' @param FDR Numeric threshold for false discovery rate to define significant markers.
#'   Default is 0.05
#' @param row_title Character string for the row title, passed to ComplexHeatmap::Heatmap.
#'   Default is ""
#' @param column_title Character string for the column title, passed to ComplexHeatmap::Heatmap.
#'   Default is ""
#' @param filterSig Logical indicating whether to filter features based on significance
#'   thresholds. Default is TRUE
#' @param propCutoff Numeric threshold for minimum proportion of cells expressing a feature.
#'   Values below this are set to 0. Default is 0.1
#' @param color Character string specifying the color for high expression values.
#'   Default is "red3"
#' @param ... Additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return A ComplexHeatmap object with bubble plot visualization
#' @export
#' @importFrom grid grid.circle
#' @importFrom ComplexHeatmap Heatmap Legend gt_render
#'
#' @examples
#' \dontrun{
#' # Get marker features
#' markerGenes <- getMarkerFeatures(
#'   ArchRProj = projHeme,
#'   useMatrix = "GeneScoreMatrix",
#'   groupBy = "Clusters"
#' )
#'
#' # Get mean and proportion matrices
#' meanMat <- getMeanMtrx(projHeme, useMatrix = "GeneScoreMatrix", groupBy = "Clusters")
#' propMat <- getNonZeroProp(projHeme, name = "GeneScoreMatrix", groupBy = "Clusters")
#'
#' # Create bubble plot
#' bubblePlot(
#'   propMtrx = propMat,
#'   meanMtrx = meanMat,
#'   markerSet = markerGenes,
#'   useFeatures = c("CD34", "CD3D", "CD8A", "CD19"),
#'   Log2FC = 1,
#'   FDR = 0.01
#' )
#' }
bubblePlot <- function(
  propMtrx = NULL,
  meanMtrx = NULL,
  markerSet = NULL,
  useFeatures = NULL,
  labelFeatures = useFeatures,
  Log2FC = 1,
  FDR = 0.05,
  row_title = "",
  column_title = "",
  filterSig = T,
  propCutoff = 0.1,
  color = "red3",
  ...
) {
  if (filterSig) {
    fdrMtrx <- assay(markerSet, "FDR") %>% as.matrix()
    rownames(fdrMtrx) <- mcols(markerSet)$name

    logfcMtrx <- assay(markerSet, "Log2FC") %>% as.matrix()
    rownames(logfcMtrx) <- mcols(markerSet)$name

    sigMtrx <- fdrMtrx < FDR & logfcMtrx > Log2FC
  }

  # if (filterSig) {
  #   useFeatures <- intersect(useFeatures, sigDf$feature)
  # }

  grps <- colnames(propMtrx)
  propMtrx <- propMtrx[useFeatures, ]
  meanMtrx <- meanMtrx[useFeatures, grps]
  sigMtrx <- sigMtrx[useFeatures, grps]

  propMtrx[propMtrx < propCutoff] <- 0

  feature2show <- ifelse(useFeatures %in% labelFeatures, useFeatures, "")
  names(feature2show) <- useFeatures

  zRange <- range(meanMtrx)
  col_fun <- circlize::colorRamp2(c(zRange[1], zRange[2]), c("white", color))

  my_layer_fun <- function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x,
      y = y,
      r = sqrt(pindex(propMtrx, i, j)) * 0.5 * min(unit.c(width, height)),
      gp = gpar(
        fill = col_fun(pindex(meanMtrx, i, j)),
        col = ifelse(pindex(sigMtrx, i, j), "black", "white")
      )
    )
  }

  message("Plotting...")

  hp <- Heatmap(
    matrix = meanMtrx,
    rect_gp = gpar(type = "none"),
    border = "black",
    col = col_fun,
    layer_fun = my_layer_fun,
    ...
  )

  lgdFill <- Legend(
    labels = as.character(seq(0.2, 1, 0.2)),
    title = "Non-zero\nProportion",
    title_gp = gpar(fontface = "plain"),
    grid_width = min(
      component_width(hp, "heatmap_body") / ncol(meanMtrx),
      component_height(hp, "heatmap_body") / nrow(meanMtrx)
    ),
    graphics = list(
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = sqrt(0.05) * min(unit.c(w, h)),
          gp = gpar(fill = "black", col = "white")
        )
      },
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = sqrt(0.1) * min(unit.c(w, h)),
          gp = gpar(fill = "black", col = "white")
        )
      },
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = sqrt(0.15) * min(unit.c(w, h)),
          gp = gpar(fill = "black", col = "white")
        )
      },
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = sqrt(0.2) * min(unit.c(w, h)),
          gp = gpar(fill = "black", col = "white")
        )
      },
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = sqrt(0.25) * min(unit.c(w, h)),
          gp = gpar(fill = "black", col = "white")
        )
      }
    )
  )

  lgdSig <- Legend(
    title = gt_render(as.character(str_glue(
      "Log<sub>2</sub>FC > {Log2FC} &<br>FDR < {FDR}"
    ))),
    at = 1,
    labels = "",
    title_gp = gpar(fontface = "plain"),
    graphics = list(
      function(x, y, w, h) {
        grid.circle(
          x,
          y,
          r = 0.4 * min(unit.c(w, h)),
          gp = gpar(col = "black", fill = "white")
        )
      }
    )
  )

  lgdList <- list(lgdFill, lgdSig)

  return(draw(
    hp,
    annotation_legend_list = lgdList,
    heatmap_legend_side = "bottom"
  ))
}
