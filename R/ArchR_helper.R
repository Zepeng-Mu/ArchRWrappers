###############################################################################
## Script purpose: Other Helper functions for analysis using ArchR
## Author: Zepeng Mu
###############################################################################

#' Add MNN adjusted reduced dimensions to an ArchR project
#'
#' This function performs Mutual Nearest Neighbors (MNN) batch correction on reduced
#' dimensions from an ArchR project and adds the corrected dimensions back to the project.
#'
#' @param ArchRProj An ArchR project object
#' @param corCutOff Numeric cutoff for correlation between library depth and reduced dimensions.
#'   Default is 0.75
#' @param reducedDims Character string specifying the name of reducedDims from ArchR project
#'   to use as input for MNN correction. Default is "IterativeLSI"
#' @param name Character string for the name of the output MNN-corrected matrix to be stored
#'   in the ArchR project. Default is "reducedMNN"
#' @param scaleDims Logical indicating whether to let ArchR scale the input reducedDims before
#'   MNN correction. Default is TRUE
#' @param scaleDimsAfter Logical indicating whether to let ArchR scale the output reducedDims
#'   after MNN correction. If NA, uses the same value as scaleDims
#' @param groupBy Character string specifying the variable name in ArchR cellColData to use
#'   as the batch variable for MNN correction
#' @param k Integer specifying the number of nearest neighbors to use in MNN correction.
#'   Default is 15
#' @param dimsToUse Numeric vector specifying which dimensions to use from the reducedDims
#'   for MNN correction
#' @param ... Additional parameters passed to \code{batchelor::reducedMNN}
#'
#' @return An ArchR project with the MNN-corrected reduced dimensions added
#' @export
#'
#' @examples
#' \dontrun{
#' projHeme <- addReducedMNN(
#'   ArchRProj = projHeme,
#'   reducedDims = "IterativeLSI",
#'   name = "reducedMNN",
#'   groupBy = "Sample",
#'   k = 15,
#'   dimsToUse = 1:30
#' )
#' }
addReducedMNN <- function(
  ArchRProj,
  corCutOff = 0.75,
  reducedDims = "IterativeLSI",
  name = "reducedMNN",
  scaleDims = T,
  scaleDimsAfter = NA,
  groupBy,
  k = 15,
  dimsToUse,
  ...
) {
  message(str_glue("Getting {reducedDims}...\n"))
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

  message("Performing reducedMNN...")
  redMnnRes <- batchelor::reducedMNN(
    origRedDim,
    k = k,
    batch = getCellColData(ArchRProj, select = groupBy, drop = T),
    ...
  )

  message(str_glue("Saving {name} to ArchR project...\n"))
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
#' This function allows you to add any custom reduced-dimensions matrix (e.g., from
#' external dimensionality reduction methods) to an ArchR project. The matrix will be
#' stored in the reducedDims slot of the project.
#'
#' @param ArchRProj An ArchR project object
#' @param name Character string specifying the name for the output reduced dimensions matrix.
#'   Default is "myReduced"
#' @param matrix A numeric matrix containing the reduced dimensions, with rows as cells
#'   (matching cell names in ArchRProj) and columns as dimensions. Must have rownames
#'   matching the cellNames in the ArchR project
#'
#' @return An ArchR project with the custom reduced dimensions matrix added
#' @export
#'
#' @examples
#' \dontrun{
#' # Add a custom UMAP embedding
#' customUMAP <- matrix(rnorm(ncol(projHeme) * 2), ncol = 2)
#' rownames(customUMAP) <- projHeme$cellNames
#' projHeme <- addAnyReducedMtrx(
#'   ArchRProj = projHeme,
#'   name = "customUMAP",
#'   matrix = customUMAP
#' )
#' }
addAnyReducedMtrx <- function(ArchRProj, name = "myReduced", matrix = NULL) {
  message(str_glue("Saving {name} to ArchR project...\n"))
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


#' Convert ArchR project to SingleCellExperiment object
#'
#' This function extracts matrices, reduced dimensions, and embeddings from an ArchR
#' project and converts them into a SingleCellExperiment object for downstream analysis
#' with Bioconductor packages.
#'
#' @param ArchRProj An ArchR project object
#' @param useMatrix Character string specifying which matrix to extract from the ArchR project.
#'   Default is "PeakMatrix"
#' @param useSeqnames Character vector of chromosome names to include in the output.
#'   Default is c("chr"%&%1:22, "chrX")
#' @param binarize Logical indicating whether to binarize the matrix (convert to presence/absence).
#'   Default is FALSE
#' @param targetAssay Character string specifying the name of the assay in the output
#'   SingleCellExperiment object. Default is "count"
#' @param scaleTo Numeric value to scale the matrix to. Default is 10000
#' @param reducedDims Character vector of reducedDims names from the ArchR project to add
#'   to the SingleCellExperiment. Default is c("IterativeLSI", "Harmony")
#' @param scaleDims Logical vector indicating whether to scale each reducedDim. Should have
#'   the same length as reducedDims. Default is c(TRUE, TRUE)
#' @param embeddings Character vector of embedding names from the ArchR project to add
#'   to the SingleCellExperiment. Default is "UMAP"
#' @param useRowData Logical indicating whether to use rowData from the ArchR matrix.
#'   If FALSE, uses rowRanges instead. Default is FALSE
#' @param threads Integer specifying the number of threads to use. Default is 4
#'
#' @return A SingleCellExperiment object containing the extracted data from the ArchR project
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' # Convert ArchR project to SingleCellExperiment
#' sce <- ArchR2sce(
#'   ArchRProj = projHeme,
#'   useMatrix = "PeakMatrix",
#'   reducedDims = c("IterativeLSI", "Harmony"),
#'   embeddings = "UMAP"
#' )
#' }
ArchR2sce <- function(
  ArchRProj,
  useMatrix = "PeakMatrix",
  useSeqnames = c("chr" %&% 1:22, "chrX"),
  binarize = F,
  targetAssay = "count",
  scaleTo = 1e4,
  reducedDims = c("IterativeLSI", "Harmony"),
  scaleDims = c(T, T),
  embeddings = "UMAP",
  useRowData = F,
  threads = 4
) {
  message(str_glue("Getting {useMatrix} from ArchR project...\n"))
  projMtrx <- getMatrixFromProject(
    ArchRProj,
    useMatrix = useMatrix,
    useSeqnames = useSeqnames,
    binarize = binarize,
    threads = threads
  )

  assayList <- list(assay(projMtrx))
  names(assayList) <- targetAssay

  rdList <- lapply(1:length(reducedDims), function(rd) {
    tmpRd <- getReducedDims(
      ArchRProj,
      reducedDims = reducedDims[rd],
      scaleDims = scaleDims[rd]
    )
    tmpRd <- tmpRd[colnames(projMtrx), ]
    return(tmpRd)
  })

  names(rdList) <- reducedDims

  embList <- lapply(1:length(embeddings), function(emb) {
    tmpEmb <- getEmbedding(ArchRProj, embedding = embeddings[emb], returnDF = T)
    tmpEmb <- tmpEmb[colnames(projMtrx), ]
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
