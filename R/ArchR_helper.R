###############################################################################
## Script purpose: Other Helper functions for analysis using ArchR
## Author: Zepeng Mu
###############################################################################

#' Add MNN adjusted reduced dimensions to an ArchR project
#'
#' @param ArchRProj An ArchR project
#' @param corCutOff Cutoff for correlation between library depth and reduced dimensions
#' @param reducedDims Name of reducedDims from ArchR project to use in reducedMNN
#' @param name Name of output MNN matrix
#' @param scaleDims Whether to let ArchR scale reducedDims
#' @param scaleDimsAfter Whether to let ArchR scale output reducedDims
#' @param groupBy Variable name in ArchR cellColData to adjust for
#' @param k 
#' @param dimsToUse Dimensions to use in reducedMNN
#' @param ... Other parameters for reducedMNN
#'
#' @return
#' @export
#'
#' @examples
addReducedMNN <- function(ArchRProj,
                          corCutOff = 0.5,
                          reducedDims = "IterativeorigRedDim",
                          name = "reducedMNN",
                          scaleDims = T,
                          scaleDimsAfter = NA,
                          groupBy,
                          k = 15,
                          dimsToUse,
                          ...) {
  cat(str_glue("Getting {reducedDims}...\n\n"))
  origRedDim <- getReducedDims(
    ArchRProj,
    corCutOff = corCutOff,
    reducedDims = reducedDims,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims
  )
  
  cat("Performing reducedMNN...\n")
  redMnnRes <- batchelor::reducedMNN(
    origRedDim,
    k = k,
    batch = getCellColData(ArchRProj, select = groupBy, drop = T),
    ...
  )
  
  cat(str_glue("Saving {name} to ArchR project...\n\n"))
  ArchRProj@reducedDims[[name]] <- SimpleList(
    matDR = redMnnRes$corrected,
    params = NA,
    date = Sys.time(),
    scaleDims = scaleDimsAfter,
    corToDepth = NA
  )
  
  return(ArchRProj)
}

#' Add any reduced-dimensions-like matrix to an ArchR project
#'
#' @param ArchRProj An ArchR project
#' @param name Name of output MNN matrix
#' @param matrix The matrix to be added to ArchR project
#'
#' @return
#' @export
#'
#' @examples
#' 
addAnyReducedMtrx <- function(ArchRProj, name = "myReduced", matrix = NULL) {
  cat(str_glue("Saving {name} to ArchR project...\n\n"))
  matrix <- matrix[ArchRProj$cellName, ]
  ArchRProj@reducedDims[[name]] <- SimpleList(
    matDR = matrix,
    params = NA,
    date = Sys.time(),
    scaleDims = NA,
    corToDepth = NA
  )
  
  return(ArchRProj)
}

#' Title
#'
#' @param ArchRProj 
#' @param matrix 
#' @param useSeqnames 
#' @param binarize 
#' @param targetAssay 
#' @param addLog 
#' @param scaleTo 
#' @param targetAssayLog 
#' @param reducedDims 
#' @param scaleDims 
#' @param threads 
#'
#' @return
#' @export
#'
#' @examples
ArchR2sce <- function(
  ArchRProj,
  matrix = "PeakMatrix",
  useSeqnames = c("chr"%&%1:22, "chrX"),
  binarize = F,
  targetAssay = "count",
  addLog = F,
  scaleTo = 1e4,
  targetAssayLog = ifelse(addLog = T, paste0("log", targetAssay)),
  reducedDims = c("IterativeLSI", "UMAP"),
  scaleDims = c(T, T),
  threads = 4
) {
  cat(str_glue("Getting {matrix} from {ArchRProj}...\n\n"))
  matrix <- getMatrixFromProject(ArchRProj, useMatrix = matrix,
                                 useSeqnames = useSeqnames,
                                 binarize = binarize, threads = threads)
  
  matrixAssay <- assay(matrix)
  if (addLog) {
    matrixAssayLog <- scale(matrixAssay, center = F,
                            scale = sparseMatrixStats::colSums2(matrixAssay))
    matrixAssayLog <- log2(matrixAssayLog + 1)
    
    assayList <- list(matrixAssay, matrixAssayLog)
    names(assayList) <- c(targetAssay, targetAssayLog)
  } else {
    assayList <- list(assay(matrix))
    names(assayList) <- targetAssay
  }
  
  rdList <- lapply(1:length(reducedDims), function(rd) {
    getReducedDims(ArchRProj, reducedDims = reducedDims[rd], scaleDims = scaleDims[rd])
  })
  
  names(rdList) <- reducedDims
  outSce <- SingleCellExperiment(
    assayList,
    colData = colData(matrix),
    rowData = rowData(matrix),
    reducedDims = rdList
  )
  
  return(outSce)
}
