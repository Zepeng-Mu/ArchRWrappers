####################################################################
# Topic Score Score Methods
####################################################################

#' Add TopicScoreMatrix to ArrowFiles or an ArchRProject
#' 
#' This function, for each sample, will independently compute sums of score for each tile
#' per cell and then infer gene topic scores.
#' This function is adapted from `addGeneScoreMatrix` function from the  `ArchR` package.
#'
#' @param genes A stranded `GRanges` object containing the ranges associated with all gene start and end coordinates. 
#' @param peaks A `GRanges` object containing the peaks with a `score` column to calculate gene-level scores. 
#' @param scoreMat 
#' @param geneModel A string giving a "gene model function" used for weighting peaks for gene score calculation. This string
#' should be a function of `x`, where `x` is the stranded distance from the transcription start site of the gene. 
#' @param matrixName The name to be used for storage of the gene activity score matrix in the provided `ArchRProject` or ArrowFiles.
#' @param extendUpstream The minimum and maximum number of basepairs upstream of the transcription start site to consider for gene
#' activity score calculation.
#' @param extendDownstream The minimum and maximum number of basepairs downstream of the transcription start site or transcription termination site 
#' (based on 'useTSS') to consider for gene activity score calculation.
#' @param useGeneBoundaries A boolean value indicating whether gene boundaries should be employed during gene activity score
#' calculation. Gene boundaries refers to the process of preventing tiles from contributing to the gene score of a given gene
#' if there is a second gene's transcription start site between the tile and the gene of interest.
#' @param geneUpstream An integer describing the number of bp upstream the gene to extend the gene body. This effectively makes the gene body larger as there
#' are proximal peaks that should be weighted equally to the gene body. This parameter is used if 'useTSS=FALSE'.
#' @param geneDownstream An integer describing the number of bp downstream the gene to extend the gene body.This effectively makes the gene body larger as there
#' are proximal peaks that should be weighted equally to the gene body. This parameter is used if 'useTSS=FALSE'.
#' @param useTSS A boolean describing whether to build gene model based on gene TSS or the gene body.
#' @param extendTSS A boolean describing whether to extend the gene TSS. By default useTSS uses the 1bp TSS while this parameter enables the extension of this
#' region with 'geneUpstream' and 'geneDownstream' respectively.
#' @param geneScaleFactor A numeric scaling factor to weight genes based on the inverse of there length i.e. [(Scale Factor)/(Gene Length)]. This
#' is scaled from 1 to the scale factor. Small genes will be the scale factor while extremely large genes will be closer to 1. This scaling helps with
#' the relative gene score value.
#' @param excludeChr A character vector containing the `seqnames` of the chromosomes that should be excluded from this analysis.
#' @param blacklist A `GRanges` object containing genomic regions to blacklist that may be extremeley over-represented and thus
#' biasing the geneScores for genes nearby that locus.
#' @param force A boolean value indicating whether to force the matrix indicated by `matrixName` to be overwritten if it already exist in the given `input`.
#' @param subThreads 
#' @param tstart 
#' @param logFile The path to a file to be used for logging ArchR output.
#'
#' @export
addTopicScoreMatrix <- function(
  genes = NULL,
  peaks = NULL,
  scoreMat = NULL,
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  peakWidth = 500,
  method = "sum",
  extendUpstream = c(1000, 100000),
  extendDownstream = c(1000, 100000),
  geneUpstream = 5000, #New Param
  geneDownstream = 0, #New Param
  useGeneBoundaries = TRUE,
  useTSS = FALSE, #New Param
  extendTSS = FALSE,
  geneScaleFactor = 5, #New Param
  excludeChr = c("chrY","chrM"),
  blacklist = NULL,
  subThreads = 1,
  tstart = NULL,
  logFile = NULL
){
  
  ArchR:::.validInput(input = genes, name = "genes", valid = c("GRanges"))
  ArchR:::.validInput(input = peaks, name = "peaks", valid = c("GRanges"))
  ArchR:::.validInput(input = scoreMat, name = "scoreMat", valid = c("matrix"))
  ArchR:::.validInput(input = geneModel, name = "geneModel", valid = c("character"))
  ArchR:::.validInput(input = extendUpstream, name = "extendUpstream", valid = c("integer"))
  ArchR:::.validInput(input = extendDownstream, name = "extendDownstream", valid = c("integer"))
  ArchR:::.validInput(input = useGeneBoundaries, name = "useGeneBoundaries", valid = c("boolean"))
  ArchR:::.validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  ArchR:::.validInput(input = blacklist, name = "blacklist", valid = c("GRanges", "null"))
  
  if(inherits(mcols(genes)$symbol, "list") | inherits(mcols(genes)$symbol, "SimpleList")){
    stop("Found a list in genes symbol! This is an incorrect format. Please correct your genes!")
  }
  
  if (is.null(peaks)) {
    stop("peaks should contain a 'name' column that can be matched to rownames in scoreMat!")
  }
  
  geneRegions <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]
  
  #Create Gene Regions Then Remove Strand Column
  if(useTSS){
    distMethod <- "GenePromoter"
    geneRegions$geneStart <- start(resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(resize(geneRegions, 1, "end"))
    geneRegions <- resize(geneRegions, 1, "start")
    if(extendTSS){
      geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)
    }
    geneRegions$geneWeight <- geneScaleFactor
  }else{
    distMethod <- "GeneBody"
    geneRegions$geneStart <- start(resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(resize(geneRegions, 1, "end"))
    geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)
    m <- 1 / width(geneRegions)
    geneRegions$geneWeight <- 1 + m * (geneScaleFactor - 1) / (max(m) - min(m))
  }
  
  #Add Gene Index For ArrowFile
  geneRegions <- sort(sortSeqlevels(geneRegions), ignore.strand = TRUE)
  
  geneRegions <- split(geneRegions, seqnames(geneRegions))
  geneRegions <- lapply(geneRegions, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })
  
  #Blacklist Split
  if(!is.null(blacklist)){
    if(length(blacklist) > 0){
      blacklist <- split(blacklist, seqnames(blacklist))
    }
  }
  
  outMat <- ArchR:::.safelapply(seq_along(geneRegions), function(z){
    chrOutMat <- tryCatch({
      #Get Gene Starts
      geneRegionz <- geneRegions[[z]]
      geneRegionz <- geneRegionz[order(geneRegionz$idx)]
      chrz <- paste0(unique(seqnames(geneRegionz)))
      
      cat(str_glue("Running {chrz}...\n\n"))
      
      # Peaks by Chr
      chrPeaks <- peaks[seqnames(peaks) == chrz]
      
      gc()
      
      # Get peak by K matrix
      chrScoreMat <- scoreMat[chrPeaks$name, ]
      
      #Time to Overlap Gene Windows
      if(useGeneBoundaries){
        
        geneStartz <- start(resize(geneRegionz, 1, "start"))
        geneEndz <- start(resize(geneRegionz, 1, "end"))
        
        pminGene <- pmin(geneStartz, geneEndz)
        pmaxGene <- pmax(geneStartz, geneEndz)
        
        idxMinus <- BiocGenerics::which(strand(geneRegionz) != "-")
        
        pReverse <- rep(max(extendDownstream), length(pminGene))
        pReverse[idxMinus] <- rep(max(extendUpstream), length(idxMinus))
        
        pReverseMin <- rep(min(extendDownstream), length(pminGene))
        pReverseMin[idxMinus] <- rep(min(extendUpstream), length(idxMinus))
        
        pForward <- rep(max(extendUpstream), length(pminGene))
        pForward[idxMinus] <- rep(max(extendDownstream), length(idxMinus))      
        
        pForwardMin <- rep(min(extendUpstream), length(pminGene))
        pForwardMin[idxMinus] <- rep(min(extendDownstream), length(idxMinus))      
        
        ################################################################
        #We will test when genes pass by another gene promoter
        ################################################################
        
        #Start of Range is based on the max observed gene ranged <- direction
        s <- pmax(
          c(1, pmaxGene[-length(pmaxGene)] + peakWidth), 
          pminGene - pReverse
        )
        s <- pmin(pminGene - pReverseMin, s)
        
        #End of Range is based on the max observed gene ranged -> direction
        e <- pmin(
          c(pminGene[-1] - peakWidth, pmaxGene[length(pmaxGene)] + pForward[length(pmaxGene)]), 
          pmaxGene + pForward
        )
        e <- pmax(pmaxGene + pForwardMin, e)
        
        extendedGeneRegion <- IRanges(start = s, end = e)
        
        idx1 <- which(pminGene - pReverseMin < start(extendedGeneRegion))
        if(length(idx1) > 0){
          stop("Error in gene boundaries minError")
        }
        
        idx2 <- which(pmaxGene + pForwardMin > end(extendedGeneRegion))
        if(length(idx2) > 0){
          stop("Error in gene boundaries maxError")
        }
        
        rm(s, e, pReverse, pReverseMin, pForward, pForwardMin, geneStartz, geneEndz, pminGene, pmaxGene)
        
      }else{
        
        extendedGeneRegion <- ranges(suppressWarnings(extendGR(geneRegionz, upstream = max(extendUpstream), downstream = max(extendDownstream))))
        
      }
      
      tmp <- suppressWarnings(IRanges::findOverlaps(extendedGeneRegion, ranges(chrPeaks)))
      x <- distance(geneRegionz[queryHits(tmp)], chrPeaks[subjectHits(tmp)])
      
      #Determine Sign for Distance relative to strand (Directionality determined based on dist from gene start)
      isMinus <- BiocGenerics::which(strand(geneRegionz) == "-")
      signDist <- sign(start(chrPeaks)[subjectHits(tmp)] - start(resize(geneRegionz,1,"start"))[queryHits(tmp)])
      signDist[isMinus] <- signDist[isMinus] * -1
      
      #Correct the orientation for the distance!
      x <- x * signDist
      
      #Evaluate Input Model
      x <- eval(parse(text = geneModel))
      
      #Get Gene Weights Related to Gene Width
      x <- x * mcols(geneRegionz)$geneWeight[queryHits(tmp)]
      
      #Remove Blacklisted Tiles!
      if(!is.null(blacklist)){
        if(length(blacklist) > 0){
          blacklistz <- blacklist[[chrz]]
          if(is.null(blacklistz) | length(blacklistz) > 0){
            peaksBlacklist <- 1 * (!overlapsAny(chrPeaks, blacklistz))
            if(sum(peaksBlacklist == 0) > 0){
              x <- x * peaksBlacklist[subjectHits(tmp)] #Multiply Such That All Blacklisted Tiles weight is now 0!
            }
          }
        }
      }
      
      #Creating Sparse Matrix
      tmp <- Matrix::sparseMatrix(
        i = queryHits(tmp),
        j = subjectHits(tmp),
        x = x,
        dims = c(length(geneRegionz), nrow(chrScoreMat))
      )
      
      if (method == "mean") {
        cntNonZero <- ncol(tmp) - sparseMatrixStats::rowCounts(tmp, 0)
        tmp <- tmp / cntNonZero
      }
      
      #Calculate Gene Scores
      chrScoreMat <- tmp %*% chrScoreMat
      rownames(chrScoreMat) <- geneRegionz$symbol
      
      #Clean Memory
      rm(isMinus, signDist, extendedGeneRegion, chrPeaks, tmp)
      gc()
      
      return(chrScoreMat)
      
    }, error = function(e){
      
      errorList <- list(
        geneRegions = geneRegions,
        blacklist = blacklist,
        chr = chrz,
        chrScoreMat = if(exists("chrScoreMat", inherits = FALSE)) chrScoreMat else "chrScoreMat"
      )
    })
    
    return(chrScoreMat)
    
  }, threads = subThreads) %>% Reduce(rbind, .)
  
  #Clean Memory
  rm(genes)
  gc()
  
  return(outMat)
  
}
