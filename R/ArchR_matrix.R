###############################################################################
## Script purpose: Helper functions for handeling ArchR matrices
## Author: Zepeng Mu
###############################################################################

#' Get mean matrix for each group from an ArchR project
#'
#' This function calculates the mean value of each feature (gene or peak) for each
#' group of cells in an ArchR project, returning a feature-by-group matrix.
#'
#' @param ArchRProj An ArchR project object
#' @param useMatrix Character string specifying which matrix to use from the ArchR project.
#'   Default is NULL (must be specified)
#' @param groupBy Character string specifying the variable in cellColData used to group cells.
#'   Default is "Sample"
#' @param scaleTo Numeric value to scale the matrix to. Default is NULL
#' @param select Character vector to select specific groups. Default is NULL (all groups)
#' @param ignoreCase Logical indicating whether to ignore case when selecting groups.
#'   Default is TRUE
#' @param threads Integer specifying the number of threads to use. Default is 4
#'
#' @return A numeric matrix with features as rows and groups as columns, containing
#'   the mean values for each feature in each group
#' @export
#'
#' @examples
#' \dontrun{
#' meanMat <- getMeanMtrx(
#'   ArchRProj = projHeme,
#'   useMatrix = "GeneScoreMatrix",
#'   groupBy = "Clusters"
#' )
#' }
getMeanMtrx <- function(ArchRProj = NULL,
                        useMatrix = NULL,
                        groupBy = "Sample",
                        scaleTo = NULL,
                        select = NULL,
                        ignoreCase = TRUE,
                        threads = 4) {
  if (!groupBy %in% colnames(getCellColData(ArchRProj))) {
    stop(stringr::str_glue("{groupBy} not in cellColData of {ArchRProj}!!"))
  }

  message("Getting matrix from ArchR project...")
  scMtrx_sce <- getMatrixFromProject(ArchRProj, useMatrix = name, useSeqnames = useSeqnames, threads = threads)
  scMtrx <- assay(scMtrx_sce)
  tmpCellCol <- getCellColData(ArchRProj, groupBy, drop = F)

  message("Calculating mean matrix...")
  meanMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpMtrx <- sparseMatrixStats::rowMeans2(scMtrx, cols = tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)

  colnames(meanMtrx) <- unique(tmpCellCol[[groupBy]])
  if (featureType == "gene") {
    rownames(meanMtrx) <- rowData(scMtrx_sce)$name
  } else if (featureType == "peak") {
    rownames(meanMtrx) <-
      stringr::str_glue("{seqnames(scMtrx_sce)}_{start(scMtrx_sce)}-{end(scMtrx_sce)}")
  }

  return(meanMtrx)
}

#' Get sum matrix for each group from an ArchR project
#'
#' This function calculates the sum of values for each feature (gene or peak) across
#' all cells in each group, returning a feature-by-group matrix.
#'
#' @param ArchRProj An ArchR project object
#' @param name Character string specifying which matrix to use from the ArchR project.
#'   Default is "PeakMatrix"
#' @param groupBy Character string specifying the variable in cellColData used to group cells.
#'   Default is "Clusters"
#' @param useSeqnames Character vector of chromosome names to include in output matrix.
#'   Default is c("chr" %&% 1:22, "chrX")
#' @param featureType Character string indicating whether matrix rows are "gene" or "peak".
#'   Default is automatically determined based on the matrix name
#' @param threads Integer specifying the number of threads to use. Default is 4
#'
#' @return A numeric matrix with features as rows and groups as columns, containing
#'   the sum of values for each feature in each group
#' @export
#'
#' @examples
#' \dontrun{
#' sumMat <- getSumMtrx(
#'   ArchRProj = projHeme,
#'   name = "PeakMatrix",
#'   groupBy = "Clusters"
#' )
#' }
getSumMtrx <- function(ArchRProj = NULL,
                       name = "PeakMatrix",
                       groupBy = "Clusters",
                       useSeqnames = c("chr" %&% 1:22, "chrX"),
                       featureType = ifelse(name == "PeakMatrix", "peak", "gene"),
                       threads = 4) {
  if (!groupBy %in% colnames(getCellColData(ArchRProj))) {
    stop(stringr::str_glue("{groupBy} not in cellColData of {ArchRProj}!!"))
  }

  message("Getting matrix from ArchR project...")
  scMtrx_sce <- getMatrixFromProject(ArchRProj, useMatrix = name, useSeqnames = useSeqnames, threads = threads)
  scMtrx <- assay(scMtrx_sce)
  tmpCellCol <- getCellColData(ArchRProj, groupBy, drop = F)

  message("Calculating sum matrix...")
  sumMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpMtrx <- sparseMatrixStats::rowSums2(scMtrx, cols = tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)

  colnames(sumMtrx) <- unique(tmpCellCol[[groupBy]])
  if (featureType == "gene") {
    rownames(sumMtrx) <- SummarizedExperiment::rowData(scMtrx_sce)$name
  } else if (featureType == "peak") {
    rownames(sumMtrx) <- SummarizedExperiment::rowData(scMtrx_sce)$name
  }

  return(sumMtrx)
}

#' Get non-zero proportion matrix from an ArchR project
#'
#' This function calculates the proportion of cells in each group that have non-zero
#' values for each feature (gene or peak), useful for identifying features with
#' sparse expression patterns.
#'
#' @param ArchRProj An ArchR project object
#' @param name Character string specifying which matrix to use from the ArchR project.
#'   Default is "GeneScoreMatrix"
#' @param groupBy Character string specifying the variable in cellColData used to group cells.
#'   Default is "Clusters"
#' @param useSeqnames Character vector of chromosome names to include in output matrix.
#'   Default is c("chr" %&% 1:22, "chrX")
#' @param threads Integer specifying the number of threads to use. Default is 14
#'
#' @return A numeric matrix with features as rows and groups as columns, containing
#'   the proportion of non-zero values (0 to 1) for each feature in each group
#' @export
#'
#' @examples
#' \dontrun{
#' propMat <- getNonZeroProp(
#'   ArchRProj = projHeme,
#'   name = "GeneScoreMatrix",
#'   groupBy = "Clusters"
#' )
#' }
getNonZeroProp <- function(ArchRProj = NULL,
                           name = "GeneScoreMatrix",
                           groupBy = "Clusters",
                           useSeqnames = c("chr" %&% 1:22, "chrX"),
                           threads = 14) {
  if (!groupBy %in% colnames(getCellColData(ArchRProj))) {
    stop(stringr::str_glue("{groupBy} not in cellColData of {ArchRProj}!!"))
  }

  message("Getting matrix from ArchR project...")
  scMtrx_sce <- getMatrixFromProject(ArchRProj, useMatrix = name, useSeqnames = useSeqnames, threads = threads)
  scMtrx <- assay(scMtrx_sce)
  tmpCellCol <- getCellColData(ArchRProj, groupBy, drop = F)

  message("Calculating non-zero proportion matrix...")
  propMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpProp <-
      sparseMatrixStats::rowCounts(scMtrx, cols = tmpCell, value = 0) / length(tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)

  propMtrx <- 1 - propMtrx
  colnames(propMtrx) <- unique(tmpCellCol[[groupBy]])
  rownames(propMtrx) <- scMtrx_sce@elementMetadata$name
  return(propMtrx)
}
