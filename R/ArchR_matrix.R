###############################################################################
## Script purpose: Helper functions for handeling ArchR matrices
## Author: Zepeng Mu
###############################################################################


#' Get mean matrix for each group from an ArchR project
#'
#' @param proj An ArchR project
#' @param name Name of matrix to use, default is GeneScoreMatrix
#' @param groupBy The variable in cellColData used to group cells, default is Clusters
#' @param useSeqnames Chromosomes to include in output matrix, default is chr1-22 and chrX
#' @param threads Number of threads to use, default is 4
#'
#' @return A feature-by-group sparse matrix
#' @export
#'
#' @examples
getMeanMtrx <- function(proj = NULL,
                        name = "GeneScoreMatrix",
                        groupBy = "Clusters",
                        useSeqnames = c("chr" %&% 1:22, "chrX"),
                        threads = 4) {
  scMtrx_summ <- getMatrixFromProject(proj, useMatrix = name, useSeqnames = useSeqnames)
  scMtrx <- assay(scMtrx_summ)
  tmpCellCol <- getCellColData(proj, groupBy, drop = F)
  meanMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpMtrx <- sparseMatrixStats::rowMeans2(scMtrx, cols = tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)
  
  colnames(meanMtrx) <- unique(tmpCellCol[[groupBy]])
  rownames(meanMtrx) <- scMtrx_summ@elementMetadata$name
  return(meanMtrx)
}

#' Get sum matrix for each group from an ArchR project
#'
#' @param proj An ArchR project
#' @param name Name of matrix to use, default is PeakMatrix
#' @param groupBy The variable in cellColData used to group cells, default is Clusters
#' @param useSeqnames Chromosomes to include in output matrix, default is chr1-22 and chrX
#' @param threads Number of threads to use, default is 4
#'
#' @return A feature-by-group sparse matrix
#' @export
#'
#' @examples
getSumMtrx <- function(proj = NULL,
                       name = "PeakMatrix",
                       groupBy = "Clusters",
                       useSeqnames = c("chr" %&% 1:22, "chrX"),
                       threads = 4) {
  scMtrx_summ <- getMatrixFromProject(proj, useMatrix = name, useSeqnames = useSeqnames)
  scMtrx <- assay(scMtrx_summ)
  tmpCellCol <- getCellColData(proj, groupBy, drop = F)
  sumMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpMtrx <- sparseMatrixStats::rowSums2(scMtrx, cols = tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)
  
  colnames(sumMtrx) <- unique(tmpCellCol[[groupBy]])
  rownames(sumMtrx) <- scMtrx_summ@elementMetadata$name
  return(sumMtrx)
}

# Get Non-zero proportion from ArchR ----
getNonZeroProp <- function(proj = NULL,
                           name = "GeneScoreMatrix",
                           groupBy = "Clusters",
                           useSeqnames = c("chr" %&% 1:22, "chrX"),
                           threads = 14) {
  scMtrx_summ <- getMatrixFromProject(proj, useMatrix = name, useSeqnames = useSeqnames)
  scMtrx <- assay(scMtrx_summ)
  tmpCellCol <- getCellColData(proj, groupBy, drop = F)
  propMtrx <- ArchR:::.safelapply(unique(tmpCellCol[[groupBy]]), function(x) {
    tmpCell <- rownames(tmpCellCol)[tmpCellCol[[groupBy]] == x]
    tmpProp <- sparseMatrixStats::rowCounts(scMtrx, cols = tmpCell, value = 0) / length(tmpCell)
  }, threads = threads) %>% Reduce("cbind", .)
  
  propMtrx <- 1 - propMtrx
  colnames(propMtrx) <- unique(tmpCellCol[[groupBy]])
  rownames(propMtrx) <- scMtrx_summ@elementMetadata$name
  return(propMtrx)
}
