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
                          reducedDims = "IterativeLSI",
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
  
  origRedDimObj <- getReducedDims(
    ArchRProj,
    corCutOff = corCutOff,
    reducedDims = reducedDims,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    returnMatrix = F
  )
  
  
  origFeatures <- origRedDimObj[[grep("Features", names(origRedDimObj))]]
  
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
    corToDepth = NA,
    Features = origFeatures,
    tileSize = 500
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
#' @param embeddings 
#' @param useRowData 
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
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
  scaleTo = 1e4,
  reducedDims = c("IterativeLSI", "Harmony"),
  scaleDims = c(T, T),
  embeddings = "UMAP",
  useRowData = F,
  threads = 4
) {
  cat(str_glue("Getting {matrix} from ArchR project...\n\n"))
  projMtrx <- getMatrixFromProject(ArchRProj, useMatrix = matrix,
                                 useSeqnames = useSeqnames,
                                 binarize = binarize, threads = threads)
  
  assayList <- list(assay(projMtrx))
  names(assayList) <- targetAssay
  
  rdList <- lapply(1:length(reducedDims), function(rd) {
    tmpRd <- getReducedDims(ArchRProj, reducedDims = reducedDims[rd], scaleDims = scaleDims[rd])
    tmpRd <- tmpRd[colnames(matrixAssay), ]
    return(tmpRd)
  })
  
  names(rdList) <- reducedDims
  
  embList <- lapply(1:length(embeddings), function(emb) {
    tmpEmb <- getEmbedding(ArchRProj, embedding = embeddings[emb], returnDF = T)
    tmpEmb <- tmpEmb[colnames(matrixAssay), ]
    return(as.matrix(tmpEmb))
  })
  
  names(embList) <- embeddings
  rdList <- append(rdList, embList)
  
  if (useRowData) {
    outSce <- SingleCellExperiment(
      assays = assayList,
      colData = colData(projMtrx),
      rowData = rowData(projMtrx),
      reducedDims = rdList
    )
  } else {
    rr <- rowRanges(projMtrx)
    names(rr) <- str_glue("{seqnames(rr)}_{start(rr)}-{end(rr)}")
    outSce <- SingleCellExperiment(
      assays = assayList,
      colData = colData(projMtrx),
      rowRanges = rr,
      reducedDims = rdList
    )
  }
  
  return(outSce)
}
