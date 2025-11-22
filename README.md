# ArchRWrappers

Utility wrappers that smooth common ArchR single-cell ATAC-seq workflows, including reduced-dimension management and group-level summaries to Seurat bridge integration and plotting helpers.

## Feature Highlights

- **Dimensionality helpers**: `addReducedMNN`, `addAnyReducedMtrx`, and `plotRedDim` automate batch correction, custom embeddings, and multi-panel visualization.
- **Cross-platform data export**: `ArchR2sce` converts ArchR projects into `SingleCellExperiment` objects for Bioconductor pipelines.
- **Matrix summaries**: `getMeanMtrx`, `getSumMtrx`, and `getNonZeroProp` compute per-group statistics that feed marker discovery or downstream plotting.
- **Visualization utilities**: `bubblePlot` and `my_theme` quickly communicate marker trends with ComplexHeatmap-backed bubbles and publication styling.
- **Integration bridge**: `addGeneBridgeIntegrationMatrix` drives Seurat/Signac bridge workflows, saving predicted RNA labels and scores back to Arrow files.
- **Convenience helpers**: `%&%`, `imessage`, and other small utilities keep scripts concise and consistent.

Documentation built from `pkgdown` lives under `docs/` (for example the reference index at `docs/reference/index.html`).

## Installation

ArchRWrappers expects a working ArchR setup (Bioconductor ≥ 3.16, R ≥ 4.2) plus the packages used inside the helpers (`ArchR`, `SingleCellExperiment`, `SummarizedExperiment`, `batchelor`, `Seurat`, `Signac`, `ComplexHeatmap`, etc.). Install the Bioconductor stack first, then pull ArchRWrappers from GitHub:

```r
# Install base dependencies once
install.packages(c("remotes", "BiocManager"))
BiocManager::install(c("ArchR", "batchelor", "SingleCellExperiment", "SummarizedExperiment"))

# Install the wrapper package
remotes::install_github("Zepeng-Mu/ArchRWrappers")
```

If you build from source, enable Arrow writing support by installing `hdf5`/`curl` system libraries as required by ArchR.

## Getting Started

### 1. Refine reduced dimensions and visualize

```r
# Batch-correct an existing embedding with batchelor and store it on the project
proj <- addReducedMNN(
	ArchRProj = proj,
	reducedDims = "IterativeLSI",
	name = "MNN",
	groupBy = "Sample"
)

# Compare embeddings side-by-side with consistent aesthetics
plotRedDim(
	ArchRProj = proj,
	reducedDims = c("IterativeLSI", "MNN"),
	colorBy = c("Sample", "Clusters"),
	savePlot = FALSE
)
```

### 2. Convert to SingleCellExperiment for Bioconductor tools

```r
sce <- ArchR2sce(
	ArchRProj = proj,
	featureSlot = "GeneScoreMatrix",
	reducedDims = c("UMAP", "MNN"),
	metadataFields = c("Sample", "Clusters")
)
```

### 3. Summarize markers with bubble plots

```r
markerGenes <- c("GATA1", "SPI1", "PAX5", "MS4A1")
meanMat <- getMeanMtrx(
	ArchRProj = proj,
	useMatrix = "GeneScoreMatrix",
	groupBy = "Clusters",
	features = markerGenes
)
propMat <- getNonZeroProp(
	ArchRProj = proj,
	useMatrix = "GeneScoreMatrix",
	groupBy = "Clusters",
	features = markerGenes
)

bubblePlot(
	meanMtx = meanMat,
	propMtx = propMat,
	name = "GeneScoreMarkers",
	rowOrder = markerGenes,
	columnOrder = colnames(meanMat)
)
```

### 4. Bridge to RNA references with Seurat

```r
proj <- addGeneBridgeIntegrationMatrix(
	ArchRProj = proj,
	matrixName = "GeneScoreMatrix",
	reducedDims = "IterativeLSI",
	seRNA = reference_rna_object,
	bridgeRNA = bridge_multiome,
	groupBy = "predictedGroup"
)
```

The call stores predicted labels (`predictedCell`, `predictedGroup`) and scores back into the ArchR project, enabling follow-up plotting or filtering without leaving ArchR.

## Contributing & Support

- Issues and pull requests are welcome on GitHub.
- See `LICENSE.md` for the MIT license terms.
- For ArchR-specific questions, check the [ArchR forum](https://github.com/GreenleafLab/ArchR) and the package documentation under `docs/`.
