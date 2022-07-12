#' Title
#'
#' @param propMtrx
#' @param meanMtrx
#' @param markerSet
#' @param useFeatures
#' @param labelFeatures
#' @param Log2FC
#' @param FDR
#' @param row_title The same as `row_title` in `Heatmap` of `ComplexHeatmap`.
#' @param column_title The same as `column_title` in `Heatmap` of `ComplexHeatmap`.
#' @param filterSig
#' @param propCutoff
#' @param color
#' @param ... Other argument passed to `Heatmap` of `ComplexHeatmap`.
#'
#' @return
#' @export
#' @importFrom grid grid.circle
#' @importFrom ComplexHeatmap Heatmap Legend gt_render
#'
#' @examples
bubblePlot <- function(propMtrx = NULL,
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
                       ...) {
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
      gp = gpar(fill = col_fun(pindex(meanMtrx, i, j)), col = ifelse(pindex(sigMtrx, i, j), "black", "white"))
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
    grid_width = min(component_width(hp, "heatmap_body") / ncol(meanMtrx), component_height(hp, "heatmap_body") / nrow(meanMtrx)),
    graphics = list(
      function(x, y, w, h) grid.circle(x, y, r = sqrt(0.05) * min(unit.c(w, h)), gp = gpar(fill = "black", col = "white")),
      function(x, y, w, h) grid.circle(x, y, r = sqrt(0.1) * min(unit.c(w, h)), gp = gpar(fill = "black", col = "white")),
      function(x, y, w, h) grid.circle(x, y, r = sqrt(0.15) * min(unit.c(w, h)), gp = gpar(fill = "black", col = "white")),
      function(x, y, w, h) grid.circle(x, y, r = sqrt(0.2) * min(unit.c(w, h)), gp = gpar(fill = "black", col = "white")),
      function(x, y, w, h) grid.circle(x, y, r = sqrt(0.25) * min(unit.c(w, h)), gp = gpar(fill = "black", col = "white"))
    )
  )

  lgdSig <- Legend(
    title = gt_render(as.character(str_glue("Log<sub>2</sub>FC > {Log2FC} &<br>FDR < {FDR}"))),
    at = 1, labels = "", title_gp = gpar(fontface = "plain"),
    graphics = list(
      function(x, y, w, h) grid.circle(x, y, r = 0.4 * min(unit.c(w, h)), gp = gpar(col = "black", fill = "white"))
    )
  )

  lgdList <- list(lgdFill, lgdSig)

  return(draw(hp, annotation_legend_list = lgdList, heatmap_legend_side = "bottom"))
}
