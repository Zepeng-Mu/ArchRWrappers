#' Calculate gene-level scores from peak-level scores
#'
#' This function calculates gene-level scores by aggregating peak-level scores based
#' on genomic distance and a gene model. It is adapted from the \code{addGeneScoreMatrix}
#' function from the ArchR package, but works with custom score matrices.
#'
#' @param genes A stranded GRanges object containing the ranges associated with all gene
#'   start and end coordinates
#' @param peaks A GRanges object containing the peaks with a 'name' column that matches
#'   rownames in scoreMat
#' @param scoreMat A numeric matrix with peaks as rows (matching peak names) and samples/topics
#'   as columns, containing the scores to aggregate
#' @param geneModel A character string giving a "gene model function" used for weighting peaks
#'   for gene score calculation. This string should be a function of \code{x}, where \code{x}
#'   is the stranded distance from the transcription start site of the gene.
#'   Default is "exp(-abs(x)/5000) + exp(-1)"
#' @param peakWidth Integer specifying the width of peaks in base pairs. Used for gene boundary
#'   calculations. Default is 500
#' @param method Character string specifying the aggregation method. Either "sum" to sum
#'   weighted peak scores, or "mean" to take the mean. Default is "sum"
#' @param extendUpstream Numeric vector of length 2 giving the minimum and maximum number of
#'   basepairs upstream of the TSS to consider for gene activity score calculation.
#'   Default is c(1000, 100000)
#' @param extendDownstream Numeric vector of length 2 giving the minimum and maximum number of
#'   basepairs downstream of the TSS (or TTS if useTSS=FALSE) to consider.
#'   Default is c(1000, 100000)
#' @param geneUpstream Integer describing the number of bp upstream the gene to extend the
#'   gene body. This effectively makes the gene body larger as there are proximal peaks
#'   that should be weighted equally to the gene body. Used if useTSS=FALSE. Default is 5000
#' @param geneDownstream Integer describing the number of bp downstream the gene to extend
#'   the gene body. Used if useTSS=FALSE. Default is 0
#' @param useGeneBoundaries Logical indicating whether gene boundaries should be employed
#'   during gene activity score calculation. Gene boundaries prevent tiles from contributing
#'   to the gene score of a given gene if there is a second gene's TSS between the tile
#'   and the gene of interest. Default is TRUE
#' @param useTSS Logical indicating whether to build gene model based on gene TSS (TRUE)
#'   or the gene body (FALSE). Default is FALSE
#' @param extendTSS Logical indicating whether to extend the gene TSS. By default useTSS
#'   uses the 1bp TSS while this parameter enables the extension of this region with
#'   geneUpstream and geneDownstream respectively. Default is FALSE
#' @param geneScaleFactor Numeric scaling factor to weight genes based on the inverse of
#'   their length, i.e., [(Scale Factor)/(Gene Length)]. This is scaled from 1 to the
#'   scale factor. Small genes will be the scale factor while extremely large genes will
#'   be closer to 1. Default is 5
#' @param excludeChr Character vector containing the seqnames of the chromosomes that should
#'   be excluded from this analysis. Default is c("chrY", "chrM")
#' @param blacklist A GRanges object containing genomic regions to blacklist that may be
#'   extremely over-represented and thus biasing the gene scores for genes nearby that locus.
#'   Default is NULL
#' @param subThreads Integer specifying the number of threads to use for parallel processing
#'   across chromosomes. Default is 1
#' @param tstart Optional POSIXct timestamp for tracking computation time. Default is NULL
#'
#' @return A numeric matrix with genes as rows and samples/topics as columns, containing
#'   the aggregated gene-level scores
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate gene scores from topic scores
#' geneScores <- scoreGeneByPeak(
#'   genes = getGenes(ArchRProj),
#'   peaks = getPeakSet(ArchRProj),
#'   scoreMat = topicScores,
#'   geneModel = "exp(-abs(x)/5000) + exp(-1)",
#'   useGeneBoundaries = TRUE
#' )
#' }
scoreGeneByPeak <- function(
  genes = NULL,
  peaks = NULL,
  scoreMat = NULL,
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  peakWidth = 500,
  method = "sum",
  extendUpstream = c(1000, 100000),
  extendDownstream = c(1000, 100000),
  geneUpstream = 5000,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  geneScaleFactor = 5,
  excludeChr = c("chrY", "chrM"),
  blacklist = NULL,
  subThreads = 1,
  tstart = NULL
) {
  ArchR:::.validInput(input = genes, name = "genes", valid = c("GRanges"))
  ArchR:::.validInput(input = peaks, name = "peaks", valid = c("GRanges"))
  ArchR:::.validInput(input = scoreMat, name = "scoreMat", valid = c("matrix"))
  ArchR:::.validInput(
    input = geneModel,
    name = "geneModel",
    valid = c("character")
  )
  ArchR:::.validInput(
    input = extendUpstream,
    name = "extendUpstream",
    valid = c("integer")
  )
  ArchR:::.validInput(
    input = extendDownstream,
    name = "extendDownstream",
    valid = c("integer")
  )
  ArchR:::.validInput(
    input = useGeneBoundaries,
    name = "useGeneBoundaries",
    valid = c("boolean")
  )
  ArchR:::.validInput(
    input = excludeChr,
    name = "excludeChr",
    valid = c("character", "null")
  )
  ArchR:::.validInput(
    input = blacklist,
    name = "blacklist",
    valid = c("GRanges", "null")
  )

  if (
    inherits(mcols(genes)$symbol, "list") |
      inherits(mcols(genes)$symbol, "SimpleList")
  ) {
    stop(
      "Found a list in genes symbol! This is an incorrect format. Please correct your genes!"
    )
  }

  if (is.null(peaks)) {
    stop(
      "peaks should contain a 'name' column that can be matched to rownames in scoreMat!"
    )
  }

  geneRegions <- genes[BiocGenerics::which(seqnames(genes) %bcni% excludeChr)]
  seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]

  # Create Gene Regions Then Remove Strand Column
  if (useTSS) {
    distMethod <- "GenePromoter"
    geneRegions$geneStart <- start(resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(resize(geneRegions, 1, "end"))
    geneRegions <- resize(geneRegions, 1, "start")
    if (extendTSS) {
      geneRegions <- extendGR(
        gr = geneRegions,
        upstream = geneUpstream,
        downstream = geneDownstream
      )
    }
    geneRegions$geneWeight <- geneScaleFactor
  } else {
    distMethod <- "GeneBody"
    geneRegions$geneStart <- start(resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(resize(geneRegions, 1, "end"))
    geneRegions <- extendGR(
      gr = geneRegions,
      upstream = geneUpstream,
      downstream = geneDownstream
    )
    m <- 1 / width(geneRegions)
    geneRegions$geneWeight <- 1 + m * (geneScaleFactor - 1) / (max(m) - min(m))
  }

  # Add Gene Index
  geneRegions <- sort(sortSeqlevels(geneRegions), ignore.strand = TRUE)

  geneRegions <- split(geneRegions, seqnames(geneRegions))
  geneRegions <- lapply(geneRegions, function(x) {
    mcols(x)$idx <- seq_along(x)
    return(x)
  })

  # Blacklist Split
  if (!is.null(blacklist)) {
    if (length(blacklist) > 0) {
      blacklist <- split(blacklist, seqnames(blacklist))
    }
  }

  # Filter peaks to those present in scoreMat
  peaksBoth <- intersect(rownames(scoreMat), peaks$name)
  peaks <- peaks[peaks$name %in% peaksBoth]
  scoreMat <- scoreMat[peaksBoth, ]

  outMat <- ArchR:::.safelapply(
    seq_along(geneRegions),
    function(z) {
      chrOutMat <- tryCatch(
        {
          # Get Gene Starts
          geneRegionz <- geneRegions[[z]]
          geneRegionz <- geneRegionz[order(geneRegionz$idx)]
          chrz <- paste0(unique(seqnames(geneRegionz)))

          message(str_glue("Running {chrz}...\n"))

          # Peaks by Chr
          chrPeaks <- peaks[seqnames(peaks) == chrz]

          gc()

          # Get peak by K matrix
          chrScoreMat <- scoreMat[chrPeaks$name, ]

          # Time to Overlap Gene Windows
          if (useGeneBoundaries) {
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
            pForwardMin[idxMinus] <- rep(
              min(extendDownstream),
              length(idxMinus)
            )

            ################################################################
            # We will test when genes pass by another gene promoter
            ################################################################

            # Start of Range is based on the max observed gene ranged <- direction
            s <- pmax(
              c(1, pmaxGene[-length(pmaxGene)] + peakWidth),
              pminGene - pReverse
            )
            s <- pmin(pminGene - pReverseMin, s)

            # End of Range is based on the max observed gene ranged -> direction
            e <- pmin(
              c(
                pminGene[-1] - peakWidth,
                pmaxGene[length(pmaxGene)] + pForward[length(pmaxGene)]
              ),
              pmaxGene + pForward
            )
            e <- pmax(pmaxGene + pForwardMin, e)

            extendedGeneRegion <- IRanges(start = s, end = e)

            idx1 <- which(pminGene - pReverseMin < start(extendedGeneRegion))
            if (length(idx1) > 0) {
              stop("Error in gene boundaries minError")
            }

            idx2 <- which(pmaxGene + pForwardMin > end(extendedGeneRegion))
            if (length(idx2) > 0) {
              stop("Error in gene boundaries maxError")
            }

            rm(
              s,
              e,
              pReverse,
              pReverseMin,
              pForward,
              pForwardMin,
              geneStartz,
              geneEndz,
              pminGene,
              pmaxGene
            )
          } else {
            extendedGeneRegion <- ranges(suppressWarnings(extendGR(
              geneRegionz,
              upstream = max(extendUpstream),
              downstream = max(extendDownstream)
            )))
          }

          tmp <- suppressWarnings(IRanges::findOverlaps(
            extendedGeneRegion,
            ranges(chrPeaks)
          ))
          x <- distance(geneRegionz[queryHits(tmp)], chrPeaks[subjectHits(tmp)])

          # Determine Sign for Distance relative to strand (Directionality determined based on dist from gene start)
          isMinus <- BiocGenerics::which(strand(geneRegionz) == "-")
          signDist <- sign(
            start(chrPeaks)[subjectHits(tmp)] -
              start(resize(geneRegionz, 1, "start"))[queryHits(tmp)]
          )
          signDist[isMinus] <- signDist[isMinus] * -1

          # Correct the orientation for the distance!
          x <- x * signDist

          # Evaluate Input Model
          x <- eval(parse(text = geneModel))

          # Get Gene Weights Related to Gene Width
          x <- x * mcols(geneRegionz)$geneWeight[queryHits(tmp)]

          # Remove Blacklisted Tiles!
          if (!is.null(blacklist)) {
            if (length(blacklist) > 0) {
              blacklistz <- blacklist[[chrz]]
              if (is.null(blacklistz) | length(blacklistz) > 0) {
                peaksBlacklist <- 1 * (!overlapsAny(chrPeaks, blacklistz))
                if (sum(peaksBlacklist == 0) > 0) {
                  x <- x * peaksBlacklist[subjectHits(tmp)] # Multiply Such That All Blacklisted Tiles weight is now 0!
                }
              }
            }
          }

          # Creating Sparse Matrix
          tmp <- Matrix::sparseMatrix(
            i = queryHits(tmp),
            j = subjectHits(tmp),
            x = x,
            dims = c(length(geneRegionz), nrow(chrScoreMat)),
            dimnames = list(geneRegionz$symbol, rownames(chrScoreMat))
          )

          # Remove empty genes and peaks
          tmp <- tmp[
            sparseMatrixStats::rowSums2(tmp) > 0,
            sparseMatrixStats::colSums2(tmp) > 0
          ]

          if (method == "mean") {
            cntNonZero <- ncol(tmp) -
              sparseMatrixStats::rowCounts(tmp, value = 0)
            tmp[cntNonZero > 0, ] <- tmp[cntNonZero > 0, ] /
              cntNonZero[cntNonZero > 0]
          }

          # Calculate Gene Scores
          chrScoreMat <- tmp %*% chrScoreMat[colnames(tmp), ]

          # Clean Memory
          rm(isMinus, signDist, extendedGeneRegion, chrPeaks, tmp)
          gc()

          return(chrScoreMat)
        },
        error = function(e) {
          errorList <- list(
            geneRegions = geneRegions,
            blacklist = blacklist,
            chr = chrz,
            chrScoreMat = if (exists("chrScoreMat", inherits = FALSE)) {
              chrScoreMat
            } else {
              "chrScoreMat"
            }
          )
        }
      )

      return(errorList)
    },
    threads = subThreads
  ) %>%
    Reduce(rbind, .)

  # Clean Memory
  rm(genes)
  gc()

  return(outMat)
}
