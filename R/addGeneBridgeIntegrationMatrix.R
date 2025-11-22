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
#'   ArchRProj = projHeme,
#'   bridge = bridgeMultiome,
#'   seRNA = seuratRNA,
#'   useMatrix = "PeakMatrix",
#'   reducedDims = "IterativeLSI",
#'   groupRNA = "celltype",
#'   addToArrow = TRUE
#' )
#'
#' # Check predicted labels
#' table(projHeme$predictedGroup)
#' }

    ArchRProj = NULL,
    bridge = NULL,
    useMatrix = "BridgePeakMatrix",
    matrixName = "GeneBridgeIntegrationMatrix",
    reducedDims = "IterativeLSI",
    annotation = NULL,
    seRNA = NULL,
    groupATAC = NULL,
    groupRNA = NULL,
    groupList = NULL,
    sampleCellsATAC = 10000,
    sampleCellsRNA = 10000,
    embeddingATAC = NULL,
    embeddingRNA = NULL,
    dimsToUse = 1:30,
    scaleDims = NULL,
    corCutOff = 0.75,
    plotUMAP = F,
    UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
    useImputation = F,
    reduction = "lsiproject",
    bridge.reduction = "cca",
    reference.reduction = "pca_harmony",
    addToArrow = F,
    scaleTo = 10000,
    genesUse = NULL,
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    transferParams = list(),
    threads = getArchRThreads(),
    verbose = TRUE,
    force = FALSE,
    logFile = createLogFile("addGeneBridgeIntegrationMatrix"),
    ...
){
  
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = matrixName, name = "matrixName", valid = c("character"))
  ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = seRNA, name = "seRNA", valid = c("SummarizedExperiment", "Seurat"))
  ArchR:::.validInput(input = groupATAC, name = "groupATAC", valid = c("character", "null"))
  ArchR:::.validInput(input = groupRNA, name = "groupRNA", valid = c("character"))
  ArchR:::.validInput(input = groupList, name = "groupList", valid = c("list", "null"))
  ArchR:::.validInput(input = sampleCellsATAC, name = "sampleCellsATAC", valid = c("integer", "null"))
  ArchR:::.validInput(input = sampleCellsRNA, name = "sampleCellsRNA", valid = c("integer", "null"))
  ArchR:::.validInput(input = embeddingATAC, name = "embeddingATAC", valid = c("data.frame", "null"))
  ArchR:::.validInput(input = embeddingRNA, name = "embeddingRNA", valid = c("data.frame", "null"))
  ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  ArchR:::.validInput(input = plotUMAP, name = "plotUMAP", valid = c("boolean"))
  ArchR:::.validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  ArchR:::.validInput(input = useImputation, name = "useImputation", valid = c("boolean"))
  ArchR:::.validInput(input = reduction, name = "reduction", valid = c("character"))
  ArchR:::.validInput(input = addToArrow, name = "addToArrow", valid = c("boolean"))
  ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  ArchR:::.validInput(input = genesUse, name = "genesUse", valid = c("character", "null"))
  ArchR:::.validInput(input = nameCell, name = "nameCell", valid = c("character"))
  ArchR:::.validInput(input = nameGroup, name = "nameGroup", valid = c("character"))
  ArchR:::.validInput(input = nameScore, name = "nameScore", valid = c("character"))
  ArchR:::.validInput(input = transferParams, name = "transferParams", valid = c("list"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  
  if(is.null(groupList)){ #If null use all cells (blocking will still occur)
    groupList <- SimpleList()
    groupList[[1]] <- SimpleList(
      ATAC = ArchRProj$cellNames,
      RNA = colnames(seRNA)
    )
  }
  
  #########################################################################################
  # 1. Check All ATAC is Accounted For!
  #########################################################################################
  message("ArchRWrappers: Checking ATAC data")
  
  if (useMatrix %ni% getAvailableMatrices(ArchRProj)) {
    stop("Matrix name provided to useMatrix does not exist in ArchRProject!")
  }
  
  if(!is.null(groupATAC)){
    dfATAC <- getCellColData(ArchRProj = ArchRProj, select = groupATAC, drop = FALSE)
  }
  nCell <- rep(0, length(ArchRProj$cellNames))
  names(nCell) <- ArchRProj$cellNames
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    ATAC <- groupList[[x]]$ATAC
    
    if(!is.null(groupATAC)){
      
      if(any(ATAC %in% dfATAC[,1])){
        idx <- which(ATAC %in% dfATAC[,1])
        ATAC2 <- rownames(dfATAC)[which(dfATAC[,1] %in% ATAC[idx])]
        if(length(idx) == length(ATAC)){
          ATAC <- ATAC2
        }else{
          ATAC <- c(ATAC[-idx], ATAC2)
        }
      }
      
    }
    
    SimpleList(ATAC = ATAC, RNA = groupList[[x]]$RNA)
    
  }) %>% SimpleList
  
  for(i in seq_along(groupList)){
    nCell[groupList[[i]]$ATAC] <- nCell[groupList[[i]]$ATAC] + 1
  }
  
  if(!all(nCell == 1)){
    stop("Missing ", length(which(nCell == 0)), " cells. Found ", length(which(nCell > 1))," overlapping cells from ArchRProj in groupList! Cannot have overlapping/missing cells in ATAC input, check 'groupList' argument!")
  }
  
  #########################################################################################
  # 2. Check All RNA is a Cell Name
  #########################################################################################
  message("ArchRWrappers: Checking seRNA data")
  
  #Set up RNA
  if(inherits(seRNA, "SummarizedExperiment")){
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if(groupRNA %ni% colnames(colData(seRNA))){
      stop("groupRNA not in colData of seRNA")
    }
    seuratRNA$Group <- paste0(colData(seRNA)[, groupRNA, drop = TRUE])
    rm(seRNA)
  }else{
    if(groupRNA %ni% colnames(seRNA@meta.data)){
      stop("groupRNA not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- paste0(seRNA@meta.data[,groupRNA])
    rm(seRNA)
  }
  
  if("RNA" %in% names(seuratRNA@assays)){
    DefaultAssay(seuratRNA) <- "RNA"
  }else{
    stop("'RNA' is not present in Seurat Object's Assays! Please make sure that this assay is present!")
  }
  gc()
  
  if(!is.null(groupRNA)){
    dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  }
  
  groupList <- lapply(seq_along(groupList), function(x){
    
    RNA <- groupList[[x]]$RNA
    
    if(!is.null(groupRNA)){
      
      if(any(RNA %in% dfRNA[,1])){
        idx <- which(RNA %in% dfRNA[,1])
        RNA2 <- rownames(dfRNA)[which(dfRNA[,1] %in% RNA[idx])]
        if(length(idx) == length(RNA)){
          RNA <- RNA2
        }else{
          RNA <- c(RNA[-idx], RNA2)
        }
      }
      
    }
    
    SimpleList(ATAC = groupList[[x]]$ATAC, RNA = RNA)
    
  }) %>% SimpleList
  
  cellRNA <- unlist(lapply(groupList, function(x) x$RNA))
  if(!all(cellRNA %in% colnames(seuratRNA))){
    stop("Found cells for RNA not in colnames(seuratRNA)! Please retry your input!")
  }
  
  seuratRNA <- seuratRNA[, unique(cellRNA)]
  
  #########################################################################################
  # 3. Create Integration Blocks
  #########################################################################################
  blockList <- SimpleList()
  
  for(i in seq_along(groupList)){
    
    gLi <- groupList[[i]]
    
    #######################################
    # ATAC
    #######################################
    message("ArchRWrappers: Making blocks for ATAC")
    
    if(length(gLi$ATAC) > sampleCellsATAC){
      
      if(!is.null(embeddingATAC)){
        probATAC <- ArchR:::.getDensity(embeddingATAC[gLi$ATAC,1], embeddingATAC[gLi$ATAC,2])$density
        probATAC <- probATAC / max(probATAC)
        cellsATAC <- gLi$ATAC[order(probATAC, decreasing = TRUE)]
      }else{
        cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      }
      
      cutoffs <- lapply(seq_len(1000), function(x) length(gLi$ATAC) / x) %>% unlist
      blockSize <- ceiling(min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1, length(gLi$ATAC)))
      
      #Density Based Blocking
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)
      
      blocks <- lapply(seq_len(nBlocks), function(x){
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }) %>% SimpleList
      
    }else{
      
      blocks <- list(gLi$ATAC)
    }
    
    #######################################
    # RNA
    #######################################
    message("ArchRWrappers: Making blocks for RNA")
    
    if(!is.null(embeddingRNA)){
      probRNA <- ArchR:::.getDensity(embeddingRNA[gLi$RNA,1], embeddingRNA[gLi$RNA,2])$density
      probRNA <- probRNA / max(probRNA)
    }else{
      probRNA <- rep(1, length(gLi$RNA))
    }
    
    blockListi <- pbmcapply::pbmclapply(seq_along(blocks), function(x){
      
      SimpleList(
        ATAC = blocks[[x]],
        RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)) , prob = probRNA)
      )
      
    }, mc.cores = 1) %>% SimpleList
    
    blockList <- c(blockList, blockListi)
    
  }
  rm(groupList)
  
  #########################################################################################
  # 4. Begin Integration
  #########################################################################################
  
  #Clean Project For Parallel
  message("ArchRWrappers: Copying ArchR Project")
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()
  
  #Gene Score Info
  peakDf <- ArchR:::.getFeatureDF(getArrowFiles(subProj), useMatrix)
  
  #Re-Index RNA
  splitPeakDf <- S4Vectors::split(peakDf, peakDf$seqnames)
  featureDF <- lapply(splitPeakDf, function(x){
    x$idx <- seq_len(nrow(x))
    return(x)
  }) %>% Reduce("rbind", .)
  dfParams <- data.frame(
    reduction = reduction
  )
  allChr <- unique(featureDF$seqnames)
  
  #Temp File Prefix
  tmpFile <- ArchR:::.tempfile()
  o <- suppressWarnings(file.remove(paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")))
  
  if(threads > 1){
    h5disableFileLocking()
  }
  
  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff,
                       dimsToUse = dimsToUse)
  
  #Create Output Directory
  message("ArchRWrappers: Creating output directories")
  
  outDir1 <- getOutputDirectory(ArchRProj)
  outDir2 <- file.path(outDir1, "RNABridgeIntegration")
  outDir3 <- file.path(outDir2, matrixName)
  dir.create(outDir1, showWarnings = FALSE)
  dir.create(outDir2, showWarnings = FALSE)
  dir.create(outDir3, showWarnings = FALSE)
  prevFiles <- list.files(outDir3, full.names = TRUE)
  prevFiles <- ArchR:::.suppressAll(file.remove(prevFiles))
  
  tstart <- Sys.time()
  
  threads2 <- max(ceiling(threads * 0.75), 1) #A Little Less here for now
  
  #Integration
  message(str_glue("ArchRWrappers: Integration in {length(blockList)} blocks\n"))
  dfAll <- ArchR:::.safelapply(seq_along(blockList), function(i){
    message(str_glue("ArchRWrappers: Integrating block {i}...\n"))
    prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
    
    blocki <- blockList[[i]]
    
    #Subset ATAC
    message("ArchRWrappers: Subsetting ATAC data")
    subProj@cellColData <- subProj@cellColData[blocki$ATAC, ]
    subProj@sampleColData <- subProj@sampleColData[unique(subProj$Sample), , drop = F]
    
    #Subset RNA
    message("ArchRWrappers: Subsetting RNA data")
    subRNA <- seuratRNA[, blocki$RNA]
    
    ##############################################################################################
    #2. Get Count Matrix, Requantify Count and Create Seurat ATAC
    ##############################################################################################
    message("ArchRWrappers: Getting partial Matrix")
    mat <- ArchR:::.getPartialMatrix(
      getArrowFiles(subProj),
      featureDF = peakDf,
      threads = 2,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = FALSE
    )
    
    rownames(mat) <- peakDf$name
    dim1 <- nrow(mat)
    mat <- mat[sparseMatrixStats::rowSums2(mat) > 0, ]
    dim2 <- nrow(mat)
    
    message(str_glue("ArchRWrappers: Removed {dim1 - dim2} peaks with 0 total counts."))
    
    atacAssay <- Signac::CreateChromatinAssay(
      counts = mat,
      sep = c(":", "-"),
      fragments = NULL,
      annotation = annotation
    )
    
    # Create Seurat sbject
    seuratATAC  <- CreateSeuratObject(counts = atacAssay, assay = "ATAC")
    
    #Set Default Assay
    DefaultAssay(seuratATAC) <- "ATAC"
    
    #Impute Matrix (its already scaled internally in ArrowFiles)
    if(useImputation){
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$reducedDims <- reducedDims
      imputeParams$dimsToUse <- dimsToUse
      imputeParams$scaleDims <- scaleDims
      imputeParams$corCutOff <- corCutOff
      imputeParams$threads <- 1
      imputeParams$logFile <- logFile
      subProj <- suppressMessages(do.call(addImputeWeights, imputeParams))
      mat <- suppressMessages(imputeMatrix(mat = mat, imputeWeights = getImputeWeights(subProj), verbose = FALSE, logFile = logFile))
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]]))) #Clean Up Space
    }
    
    # Normalize query
    seuratATAC <- Signac::RunTFIDF(seuratATAC, assay = "ATAC")
    
    #Clean Memory
    rm(mat)
    gc()
    
    # Drop first dimension for ATAC reduction
    dims.atac <- 2:50
    dims.rna <- 1:50
    
    subRNA <- SCTransform(object = subRNA, conserve.memory = T, ncells = 3000) %>%
      RunPCA() %>%
      RunUMAP(dims = 1:50, return.model = T)
    
    DefaultAssay(subRNA) <- "SCT"
    
    subRNAExt <- PrepareBridgeReference(
      reference = subRNA,
      bridge = bridge,
      reference.reduction = "pca",
      reference.dims = 1:50,
      normalization.method = "SCT",
      verbose = T
    )
    
    ##############################################################################################
    #3. Transfer Anchors
    ##############################################################################################
    message("ArchRWrappers: Transfering Anchors")
    transferAnchors <- Seurat::FindBridgeTransferAnchors(
      extended.reference = subRNAExt,
      query = seuratATAC,
      reduction = reduction,
      bridge.reduction = bridge.reduction,
      dims = dims.atac
    )
    
    ##############################################################################################
    #4. Transfer Data
    ##############################################################################################
    rDSub <- rD[colnames(seuratATAC), , drop = F]
    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = rDSub,
      key = "LSI_",
      assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- F
    transferParams$dims <- seq_len(ncol(rDSub))
    
    #Group
    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)
    
    #RNA Names
    # transferParams$refdata <- colnames(subRNA)
    # rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[, 1]
    rnaLabels2 <- rownames(rnaLabels)
    
    if(addToArrow){
      transferParams$refdata <- GetAssayData(subRNA, assay = "RNA", slot = "data")
      gc()
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }
    
    #Match results
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC),
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id,
      predictedCell = rnaLabels2
    )
    rownames(matchDF) <- matchDF$cellNames
    
    # jointCCA <- DataFrame(transferAnchors@object.list[[1]]@reductions$cca@cell.embeddings)
    # jointCCA$Assay <- ifelse(endsWith(rownames(jointCCA), "_reference"), "RNA", "ATAC")
    # jointCCA$Group <- NA
    # jointCCA$Score <- NA
    # jointCCA[paste0(colnames(subRNA), "_reference"), "Group"] <- subRNA$Group
    # jointCCA[paste0(matchDF$cellNames, "_query"), "Group"] <- matchDF$predictedGroup
    # jointCCA[paste0(matchDF$cellNames, "_query"), "Score"] <- matchDF$predictionScore
    # ArchR:::.safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
    
    #Clean Memory
    rm(transferParams, transferAnchors)
    gc()
    
    ##############################################################################################
    #5. Add To Temp Hdf5
    ##############################################################################################
    
    if(addToArrow){
      
      #Quickly Write to A Temp Hdf5 File Split By Sample to Then Enable Writing to Each Arrow File
      
      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- .safeSubset( #If Rownames disappeared this will catch that!
        mat = matchedRNA,
        subsetRows = paste0(featureDF$name),
        subsetCols = matchDF$cellNames
      )
      
      for(z in seq_along(uniqueSamples)){
        
        mat <- matchedRNA[, which(sampleNames == uniqueSamples[z]), drop = FALSE]
        Group <- uniqueSamples[z]
        
        o <- tryCatch({h5delete(tmpFilei, paste0(Group))}, error = function(x){})
        o <- h5createGroup(tmpFilei, paste0(Group))
        
        #Convert Columns to Rle
        j <- Rle(findInterval(seq(mat@x)-1, mat@p[-1]) + 1)
        
        #Info
        lengthRle <- length(j@lengths)
        lengthI <- length(mat@i)
        
        #Create Data Set
        o <- ArchR:::.suppressAll(h5createDataset(tmpFilei, paste0(Group,"/i"), storage.mode = "integer",
                                                  dims = c(lengthI, 1), level = 0))
        
        o <- ArchR:::.suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jLengths"), storage.mode = "integer",
                                                  dims = c(lengthRle, 1), level = 0))
        
        o <- ArchR:::.suppressAll(h5createDataset(tmpFilei, paste0(Group,"/jValues"), storage.mode = "integer",
                                                  dims = c(lengthRle, 1), level = 0))
        
        o <- ArchR:::.suppressAll(h5createDataset(tmpFilei, paste0(Group, "/x"), storage.mode = "double",
                                                  dims = c(lengthI, 1), level = 0))
        
        #Write Data Set
        o <- ArchR:::.suppressAll(h5write(obj = mat@i + 1, file = tmpFilei, name = paste0(Group,"/i")))
        o <- ArchR:::.suppressAll(h5write(obj = j@lengths, file = tmpFilei, name = paste0(Group,"/jLengths")))
        o <- ArchR:::.suppressAll(h5write(obj = j@values, file = tmpFilei, name = paste0(Group,"/jValues")))
        o <- ArchR:::.suppressAll(h5write(obj = mat@x, file = tmpFilei, name = paste0(Group, "/x")))
        o <- ArchR:::.suppressAll(h5write(obj = colnames(mat), file = tmpFilei, name = paste0(Group, "/cellNames")))
        #Row Names is always the same
        
      }
      
      rm(matchedRNA, mat, j)
      
    }
    
    gc()
    
    matchDF$Block <- Rle(i)
    return(matchDF)
    
  }, threads = threads2) %>% Reduce("rbind", .)
  
  ##############################################################################################
  #5. Plot UMAPs for Co-Embeddings from CCA
  ##############################################################################################
  if(plotUMAP){
    
    for(i in seq_along(blockList)){
      
      o <- tryCatch({
        
        prefix <- sprintf("Block (%s of %s) :", i , length(blockList))
        
        jointCCA <- readRDS(file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        set.seed(1) # Always do this prior to UMAP
        UMAPParams <- .mergeParams(UMAPParams, list(n_neighbors = 40, min_dist = 0.4, metric="cosine", verbose=FALSE))
        UMAPParams$X <- as.data.frame(jointCCA[, grep("CC_", colnames(jointCCA))])
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- tryCatch({
          do.call(uwot::umap, UMAPParams)
        }, error = function(e){
          errorList <- UMAPParams
          ArchR:::.logError(e, fn = "uwot::umap", info = prefix, errorList = errorList, logFile = logFile)
        })
        
        #Add UMAP and Save Again
        jointCCA$UMAP1 <- uwotUmap[,1]
        jointCCA$UMAP2 <- uwotUmap[,2]
        .safeSaveRDS(object = jointCCA, file = file.path(outDir3, paste0("Save-Block", i,"-JointCCA.rds")))
        
        p1 <- ggPoint(
          x = uwotUmap[,1],
          y = uwotUmap[,2],
          color = jointCCA$Assay,
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        p2 <- ggPoint(
          x = uwotUmap[,1],
          y = uwotUmap[,2],
          color = jointCCA$Group,
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by scRNA Group"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        pdf(file.path(outDir3, paste0("Save-Block", i,"-JointCCA-UMAP.pdf")), width = 12, height = 6, useDingbats = FALSE)
        ggAlignPlots(p1,p2,type="h")
        dev.off()
        
      }, error = function(e){
        
      })
      
    }
    
  }
  
  ##############################################################################################
  #6. Read sub-matrices and store in ArrowFiles
  ##############################################################################################
  
  if(addToArrow){
    
    matrixName <- ArchR:::.isProtectedArray(matrixName)
    
    integrationFiles <- paste0(tmpFile, "-IntegrationBlock-", seq_along(blockList), ".h5")
    
    if(!all(file.exists(integrationFiles))){
      stop("Something went wrong with integration as not all temporary files containing integrated RNA exist!")
    }
    
    h5list <- ArchR:::.safelapply(seq_along(integrationFiles), function(x){
      h5ls(integrationFiles[x])
    }, threads = threads)
    
    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)
    
    o <- ArchR:::.safelapply(seq_along(allSamples), function(y){
      
      sample <- allSamples[y]
      
      prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))
      
      sampleIF <- lapply(seq_along(h5list), function(x){
        if(any(h5list[[x]]$group==paste0("/",sample))){
          integrationFiles[x]
        }else{
          NULL
        }
      }) %>% unlist
      
      sampleMat <- lapply(seq_along(sampleIF), function(x){
        
        cellNames <- .h5read(sampleIF[x], paste0(sample, "/cellNames"))
        
        mat <- sparseMatrix(
          i = ArchR:::.h5read(sampleIF[x], paste0(sample, "/i"))[,1],
          j = as.vector(
            Rle(
              ArchR:::.h5read(sampleIF[x], paste0(sample, "/jValues"))[,1],
              ArchR:::.h5read(sampleIF[x], paste0(sample, "/jLengths"))[,1]
            )
          ),
          x = ArchR:::.h5read(sampleIF[x], paste0(sample, "/x"))[,1],
          dims = c(nrow(featureDF), length(cellNames))
        )
        colnames(mat) <- cellNames
        
        mat
        
      }) %>% Reduce("cbind", .)
      
      sampleMat@x <- exp(sampleMat@x) - 1 #Back To Counts
      sampleMat <- ArchR:::.normalizeCols(sampleMat, scaleTo = scaleTo) #Scale to 10,000
      sampleMat <- drop0(sampleMat) # Drop 0's
      rownames(sampleMat) <- paste0(featureDF$name)
      sampleMat <- sampleMat[,ArchRProj$cellNames[BiocGenerics::which(ArchRProj$Sample == sample)], drop = FALSE]
      
      ######################################
      # Initialize SP Mat Group
      ######################################
      o <- ArchR:::.createArrowGroup(ArrowFile = ArrowFiles[sample], group = matrixName, force = force)
      
      o <- ArchR:::.initializeMat(
        ArrowFile = ArrowFiles[sample],
        Group = matrixName,
        Class = "double",
        Units = "NormCounts",
        cellNames = colnames(sampleMat),
        params = dfParams,
        featureDF = featureDF,
        force = force
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictionScore"],
        file = ArrowFiles[sample],
        name = paste0(matrixName, "/Info/predictionScore")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedGroup"],
        file = ArrowFiles[sample],
        name = paste0(matrixName, "/Info/predictedGroup")
      )
      
      o <- h5write(
        obj = dfAll[colnames(sampleMat), "predictedCell"],
        file = ArrowFiles[sample],
        name = paste0(matrixName, "/Info/predictedCell")
      )
      
      for(z in seq_along(allChr)){
        
        chrz <- allChr[z]
        
        idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
        matz <- sampleMat[idz, ,drop=FALSE]
        stopifnot(identical(paste0(featureDF$name[idz]), paste0(rownames(matz))))
        
        #Write sparseMatrix to Arrow File!
        o <- ArchR:::.addMatToArrow(
          mat = matz,
          ArrowFile = ArrowFiles[sample],
          Group = paste0(matrixName, "/", chrz),
          binarize = FALSE,
          addColSums = TRUE,
          addRowSums = TRUE,
          addRowVarsLog2 = TRUE,
          logFile = logFile
        )
        
        #Clean Memory
        rm(matz)
        
        if(z %% 3 == 0 | z == length(allChr)){
          gc()
        }
        
      }
      
      0
      
    }, threads = threads)
    
    o <- suppressWarnings(file.remove(integrationFiles))
    
  }
  
  message("ArchRWrappers: Adding cellColData")
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictedCell,
    name = nameCell,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictedGroup,
    name = nameGroup,
    force = TRUE
  )
  
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictionScore,
    name = nameScore,
    force = TRUE
  )
  
  return(ArchRProj)
  
}
