#' Integrate ATAC-seq and RNA-seq data using Seurat Bridge Integration
#'
#' This function performs bridge integration between an ArchR project (ATAC-seq) and
#' a Seurat object or SummarizedExperiment (RNA-seq) using a bridge dataset. It uses
#' Seurat's bridge integration workflow to transfer labels and optionally imputed gene
#' expression from RNA to ATAC cells.
#'
#' @param ArchRProj An ArchR project object containing ATAC-seq data
#' @param bridge A Seurat object containing the bridge dataset (e.g., multiome data)
#' @param useMatrix Character string specifying which matrix to use from the ArchR project.
#'   Default is "BridgePeakMatrix"
#' @param matrixName Character string for the name of the output integration matrix.
#'   Default is "GeneBridgeIntegrationMatrix"
#' @param reducedDims Character string specifying which reducedDims to use from the ArchR
#'   project. Default is "IterativeLSI"
#' @param annotation A GRanges object with gene annotations for creating ChromatinAssay.
#'   Default is NULL
#' @param seRNA A Seurat object or SummarizedExperiment containing RNA-seq data to integrate
#' @param groupATAC Character string specifying the grouping variable in ArchR cellColData
#'   for ATAC cells. Default is NULL
#' @param groupRNA Character string specifying the grouping variable in seRNA metadata
#'   for RNA cells
#' @param groupList A list of SimpleList objects defining cell groupings for integration.
#'   Each element should have ATAC and RNA components. Default is NULL (use all cells)
#' @param sampleCellsATAC Integer specifying the maximum number of ATAC cells to sample
#'   per block. Default is 10000
#' @param sampleCellsRNA Integer specifying the maximum number of RNA cells to sample
#'   per block. Default is 10000
#' @param embeddingATAC A data.frame with ATAC cell embeddings for density-based sampling.
#'   Default is NULL
#' @param embeddingRNA A data.frame with RNA cell embeddings for density-based sampling.
#'   Default is NULL
#' @param dimsToUse Numeric vector specifying which dimensions to use from reducedDims.
#'   Default is 1:30
#' @param scaleDims Logical indicating whether to scale dimensions. Default is NULL
#' @param corCutOff Numeric cutoff for correlation filtering. Default is 0.75
#' @param plotUMAP Logical indicating whether to plot UMAP visualizations of the integration.
#'   Default is FALSE
#' @param UMAPParams List of parameters for UMAP calculation. Default is
#'   list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE)
#' @param useImputation Logical indicating whether to use ArchR's imputation weights.
#'   Default is FALSE
#' @param reduction Character string specifying the reduction method for ATAC data.
#'   Default is "lsiproject"
#' @param bridge.reduction Character string specifying the bridge reduction method.
#'   Default is "cca"
#' @param reference.reduction Character string specifying the reference reduction method.
#'   Default is "pca_harmony"
#' @param addToArrow Logical indicating whether to add the integrated gene expression
#'   matrix to Arrow files. Default is FALSE
#' @param scaleTo Numeric value to scale the integrated matrix to. Default is 10000
#' @param genesUse Character vector of genes to use for integration. Default is NULL (all genes)
#' @param nameCell Character string for the name of the predicted cell column in cellColData.
#'   Default is "predictedCell"
#' @param nameGroup Character string for the name of the predicted group column in cellColData.
#'   Default is "predictedGroup"
#' @param nameScore Character string for the name of the prediction score column in cellColData.
#'   Default is "predictedScore"
#' @param transferParams List of additional parameters to pass to Seurat::TransferData.
#'   Default is list()
#' @param threads Integer specifying the number of threads to use. Default is getArchRThreads()
#' @param verbose Logical indicating whether to print verbose messages. Default is TRUE
#' @param force Logical indicating whether to force overwrite of existing data. Default is FALSE
#' @param logFile Character string specifying the path to the log file.
#'   Default is createLogFile("addGeneBridgeIntegrationMatrix")
#' @param ... Additional arguments (currently unused)
#'
#' @return An ArchR project with predicted cell labels, groups, and scores added to cellColData.
#'   If addToArrow=TRUE, also adds the integrated gene expression matrix to Arrow files.
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform bridge integration
#' projHeme <- addGeneBridgeIntegrationMatrix(
#'     ArchRProj = projHeme,
#'     bridge = bridgeMultiome,
#'     seRNA = seuratRNA,
#'     useMatrix = "PeakMatrix",
#'     reducedDims = "IterativeLSI",
#'     groupRNA = "celltype",
#'     addToArrow = TRUE
#' )
#'
#' # Check predicted labels
#' table(projHeme$predictedGroup)
#' }
