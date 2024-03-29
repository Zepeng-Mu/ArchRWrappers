###############################################################################
## Script purpose: Helper functions for handeling ArchR matrices
## Author: Zepeng Mu
###############################################################################

#' Get mean matrix for each group from an ArchR project
#'
#' @param ArchRProj An ArchR project
#' @param name Name of matrix to use, default is GeneScoreMatrix
#' @param groupBy The variable in cellColData used to group cells, default is Clusters
#' @param useSeqnames Chromosomes to include in output matrix, default is chr1-22 and chrX
#' @param threads Number of threads to use, default is 4
#'
#' @return A feature-by-group sparse matrix
#' @export
#'
#' @examples
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
#' @param ArchRProj An ArchR project
#' @param name Name of matrix to use, default is PeakMatrix
#' @param groupBy The variable in cellColData used to group cells, default is Clusters
#' @param useSeqnames Chromosomes to include in output matrix, default is chr1-22 and chrX
#' @param threads Number of threads to use, default is 4
#' @param featureType String indicates whether matrix rows are `gene` or `peak`, or anything else
#'
#' @return A feature-by-group sparse matrix
#' @export
#'
#' @examples
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

#' Get Non-zero proportion from ArchR project
#'
#' @param ArchRProj An ArchR project
#' @param name Name of matrix to use, default is GeneScoreMatrix
#' @param groupBy The variable in cellColData used to group cells, default is Clusters
#' @param useSeqnames Chromosomes to include in output matrix, default is chr1-22 and chrX
#' @param threads Number of threads to use, default is 4
#'
#' @return
#' @export
#'
#' @examples
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
