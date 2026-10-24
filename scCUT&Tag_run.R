library(chromVAR)
library(tidyverse)
library(ArchR)
library(parallel)
library(Seurat)

set.seed(1)
addArchRThreads(threads = 8)
addArchRGenome("hg38")


inputFiles_k4me3 <- c(
  "sck4me3-w6-e1-fragmentsInCells.tsv.gz",
  "sck4me3-w6-e2-fragmentsInCells.tsv.gz",
  "sck4me3-w6-e3-fragmentsInCells.tsv.gz"
)
names(inputFiles_k4me3) <- c("humanw6-k4me3-1", "humanw6-k4me3-2", "humanw6-k4me3-3")
ArrowFiles_k4me3 <- createArrowFiles(
  inputFiles = inputFiles_k4me3,
  sampleNames = names(inputFiles_k4me3),
  minTSS = 0,
  minFrags = 200,
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 5000),
  promoterRegion = c(5000, 0),
  TSSParams = list(flank = 5000),
  force = FALSE,
  addGeneScoreMat = TRUE,
  threads = 8
)

doubScoresk27_k4me3 <- addDoubletScores(
  input = ArrowFiles_k4me3,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 8,
  force = FALSE
)

scK4me3w6.new.proj <- ArchRProject(
  ArrowFiles = ArrowFiles_k4me3,
  outputDirectory = "w6k4me3",
  copyArrows = TRUE
)

scK4me3w6.new.proj$`log10(nFrags)` <- log10(scK4me3w6.new.proj$nFrags)
getAvailableMatrices(scK4me3w6.new.proj)
scK4me3w6.new.proj
scK4me3w6.new.proj.rmdoub <- filterDoublets(ArchRProj = scK4me3w6.new.proj)
scK4me3w6.new.proj.rmdoub

plotGroups(
  scK4me3w6.new.proj.rmdoub,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) +
  labs(x = NULL, y = "nFrags") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(400, 1.5)
  )
plotFragmentSizes(ArchRProj = scK4me3w6.new.proj.rmdoub, maxSize = 2000) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(.7, .7),
    legend.text = element_text(size = 10)
  )
plotTSSEnrichment(ArchRProj = scK4me3w6.new.proj.rmdoub) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.text = element_text(size = 10)
  )

scK4me3w6.new.proj.rmdoub <- addIterativeLSI(
  ArchRProj = scK4me3w6.new.proj.rmdoub,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  threads = 8,
  iterations = 2,
  clusterParams = list(
    resolution = 0.5,
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE,
  LSIMethod = 2,
  scaleDims = TRUE,
  nPlot = 10000,
  varFeatures = 25000,
  selectionMethod = "var",
  dimsToUse = 1:30
)
scK4me3w6.new.proj.rmdoub <- addHarmony(
  ArchRProj = scK4me3w6.new.proj.rmdoub,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  force = TRUE,
  groupBy = "Sample"
)
scK4me3w6.new.proj.rmdoub <- addClusters(
  input = scK4me3w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  force = TRUE,
  resolution = 0.5,
  maxClusters = 100,
  dimsToUse = 2:30
)
scK4me3w6.new.proj.rmdoub <- addUMAP(
  ArchRProj = scK4me3w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  name = "IterativeLSIUMAP",
  force = TRUE,
  threads = 8,
  dimsToUse = 2:30,
  nNeighbors = 40,
  minDist = 0.7
)
p2ac <- plotEmbedding(
  ArchRProj = scK4me3w6.new.proj.rmdoub,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "IterativeLSIUMAP",
  log2Norm = FALSE,
  rastr = FALSE,
  size = 0.5,
  labelMeans = TRUE
)
p2ac

scK4me3w6.new.proj.rmdoub <- addImputeWeights(scK4me3w6.new.proj.rmdoub)

table(scK4me3w6.new.proj.rmdoub$Clusters, useNA = "ifany")
scK4me3w6.new.proj.rmdoub.markersGS <- getMarkerFeatures(
  ArchRProj = scK4me3w6.new.proj.rmdoub,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  k = 100,
  maxCells = 500,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads = 8
)

plotEmbedding(
  scK4me3w6.new.proj.rmdoub,
  name = "cellname",
  colorBy = "cellColData",
  pal = pals,
  embedding = "IterativeLSIUMAP",
  labelMeans = TRUE
) +
  labs(title = NULL)
plotEmbedding(
  scK4me3w6.new.proj.rmdoub,
  name = "Sample",
  colorBy = "cellColData",
  rastr = TRUE,
  size = 0.15,
  randomize = TRUE,
  seed = 404,
  bins = 1000,
  embedding = "IterativeLSIUMAP",
  labelMeans = FALSE
)


inputFiles <- c(
  "./sck27ac-w6-e1-fragmentsInCells.tsv.gz",
  "./sck27ac-w6-e2-fragmentsInCells.tsv.gz",
  "./sck27ac-w6-e3-fragmentsInCells.tsv.gz"
)
names(inputFiles) <- c("humanw6-1", "humanw6-2", "humanw6-3")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 0,
  minFrags = 200,
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 5000),
  promoterRegion = c(5000, 0),
  TSSParams = list(flank = 5000),
  force = TRUE,
  addGeneScoreMat = TRUE,
  threads = 8
)

doubScoresk27 <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 8,
  force = TRUE
)

scK27w6.1.proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "humanw8",
  copyArrows = TRUE
)

scK27w6.1.proj$`log10(nFrags)` <- log10(scK27w6.1.proj$nFrags)
getAvailableMatrices(scK27w6.1.proj)
scK27w6.1.proj
scK27w6.1.proj.rmdoub <- filterDoublets(ArchRProj = scK27w6.1.proj)
scK27w6.1.proj.rmdoub

plotFragmentSizes(ArchRProj = scK27w6.1.proj.rmdoub)
plotTSSEnrichment(ArchRProj = scK27w6.1.proj.rmdoub)

scK27w6.1.proj.rmdoub <- addIterativeLSI(
  ArchRProj = scK27w6.1.proj.rmdoub,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  threads = 8,
  iterations = 2,
  clusterParams = list(
    resolution = 0.5,
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE,
  LSIMethod = 2,
  scaleDims = FALSE,
  nPlot = 10000,
  varFeatures = 50000,
  selectionMethod = "var",
  dimsToUse = 1:30
)
scK27w6.1.proj.rmdoub <- addHarmony(
  ArchRProj = scK27w6.1.proj.rmdoub,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  force = TRUE,
  groupBy = "Sample"
)
scK27w6.1.proj.rmdoub <- addClusters(
  input = scK27w6.1.proj.rmdoub,
  reducedDims = "Harmony",
  force = TRUE,
  resolution = 0.5,
  maxClusters = 100,
  dimsToUse = 2:30
)
scK27w6.1.proj.rmdoub <- addUMAP(
  ArchRProj = scK27w6.1.proj.rmdoub,
  reducedDims = "Harmony",
  name = "IterativeLSIUMAP",
  force = TRUE,
  threads = 8,
  dimsToUse = 2:30,
  nNeighbors = 40,
  minDist = 0.7
)

p2ac <- plotEmbedding(
  ArchRProj = scK27w6.1.proj.rmdoub,
  colorBy = "cellColData",
  name = "cellname",
  embedding = "IterativeLSIUMAP",
  log2Norm = FALSE,
  rastr = FALSE,
  size = 0.5,
  labelMeans = TRUE
)
p2ac
p2ac <- plotEmbedding(
  ArchRProj = scK27w6.1.proj.rmdoub,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "IterativeLSIUMAP",
  log2Norm = FALSE,
  rastr = FALSE,
  size = 0.5,
  labelMeans = TRUE
)
p2ac

inputFiles_k27me3 <- c(
  "sck27me3-w6-e1-fragmentsInCells.tsv.gz",
  "sck27me3-w6-e2-fragmentsInCells.tsv.gz",
  "sck27me3-w6-e3-fragmentsInCells.tsv.gz"
)
names(inputFiles_k27me3) <- c("humanw6-k27me3-1", "humanw6-k27me3-2", "humanw6-k27me3-3")
ArrowFiles_k27me3 <- createArrowFiles(
  inputFiles = inputFiles_k27me3,
  sampleNames = names(inputFiles_k27me3),
  minTSS = 0,
  minFrags = 200,
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 5000),
  promoterRegion = c(5000, 0),
  TSSParams = list(flank = 5000),
  force = FALSE,
  addGeneScoreMat = TRUE,
  threads = 8
)

doubScoresk27_k27me3 <- addDoubletScores(
  input = ArrowFiles_k27me3,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 8,
  force = FALSE
)

scK27me3w6.new.proj <- ArchRProject(
  ArrowFiles = ArrowFiles_k27me3,
  outputDirectory = "w6k27me3",
  copyArrows = TRUE
)

scK27me3w6.new.proj$`log10(nFrags)` <- log10(scK27me3w6.new.proj$nFrags)
getAvailableMatrices(scK27me3w6.new.proj)
scK27me3w6.new.proj
scK27me3w6.new.proj.rmdoub <- filterDoublets(ArchRProj = scK27me3w6.new.proj)
scK27me3w6.new.proj.rmdoub

plotGroups(
  scK27me3w6.new.proj.rmdoub,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) +
  labs(x = NULL, y = "nFrags") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(400, 1.5)
  )
plotFragmentSizes(ArchRProj = scK27me3w6.new.proj.rmdoub, maxSize = 2000) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(.7, .7),
    legend.text = element_text(size = 10)
  )

scK27me3w6.new.proj.rmdoub <- addIterativeLSI(
  ArchRProj = scK27me3w6.new.proj.rmdoub,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  threads = 8,
  iterations = 2,
  clusterParams = list(
    resolution = 0.5,
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE,
  LSIMethod = 2,
  scaleDims = TRUE,
  nPlot = 10000,
  varFeatures = 25000,
  selectionMethod = "var",
  dimsToUse = 1:30
)
scK27me3w6.new.proj.rmdoub <- addHarmony(
  ArchRProj = scK27me3w6.new.proj.rmdoub,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  force = TRUE,
  groupBy = "Sample"
)
scK27me3w6.new.proj.rmdoub <- addClusters(
  input = scK27me3w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  force = TRUE,
  resolution = 1,
  maxClusters = 100,
  dimsToUse = 2:30
)
scK27me3w6.new.proj.rmdoub <- addUMAP(
  ArchRProj = scK27me3w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  name = "IterativeLSIUMAP",
  force = TRUE,
  threads = 8,
  dimsToUse = 2:30,
  nNeighbors = 40,
  minDist = 0.7
)
p2ac <- plotEmbedding(
  ArchRProj = scK27me3w6.new.proj.rmdoub,
  colorBy = "cellColData",
  name = "cellname",
  embedding = "IterativeLSIUMAP",
  log2Norm = FALSE,
  rastr = FALSE,
  size = 0.5,
  labelMeans = TRUE
)
p2ac

inputFiles_k4me1 <- c(
  "sck4me1-w6-e1-fragmentsInCells.tsv.gz",
  "sck4me1-w6-e2-fragmentsInCells.tsv.gz",
  "sck4me1-w6-e3-fragmentsInCells.tsv.gz"
)
names(inputFiles_k4me1) <- c("humanw6-k4me1-1", "humanw6-k4me1-2", "humanw6-k4me1-3")
ArrowFiles_k4me1 <- createArrowFiles(
  inputFiles = inputFiles_k4me1,
  sampleNames = names(inputFiles_k4me1),
  minTSS = 0,
  minFrags = 200,
  addTileMat = TRUE,
  TileMatParams = list(tileSize = 5000),
  promoterRegion = c(5000, 0),
  TSSParams = list(flank = 5000),
  force = FALSE,
  addGeneScoreMat = TRUE,
  threads = 8
)

doubScoresk27_k4me1 <- addDoubletScores(
  input = ArrowFiles_k4me1,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  threads = 8,
  force = FALSE
)

scK4me1w6.new.proj <- ArchRProject(
  ArrowFiles = ArrowFiles_k4me1,
  outputDirectory = "w6k4me1",
  copyArrows = TRUE
)

scK4me1w6.new.proj$`log10(nFrags)` <- log10(scK4me1w6.new.proj$nFrags)
getAvailableMatrices(scK4me1w6.new.proj)
scK4me1w6.new.proj
scK4me1w6.new.proj.rmdoub <- filterDoublets(ArchRProj = scK4me1w6.new.proj)
scK4me1w6.new.proj.rmdoub


plotGroups(
  scK4me1w6.new.proj.rmdoub,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) +
  labs(x = NULL, y = "nFrags") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(400, 1.5)
  )
plotFragmentSizes(ArchRProj = scK4me1w6.new.proj.rmdoub, maxSize = 2000) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = c(.7, .7),
    legend.text = element_text(size = 10)
  )

scK4me1w6.new.proj.rmdoub <- addIterativeLSI(
  ArchRProj = scK4me1w6.new.proj.rmdoub,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  threads = 8,
  iterations = 2,
  clusterParams = list(
    resolution = 0.5,
    sampleCells = 10000,
    n.start = 10
  ),
  force = TRUE,
  LSIMethod = 2,
  scaleDims = TRUE,
  nPlot = 10000,
  varFeatures = 25000,
  selectionMethod = "var",
  dimsToUse = 1:30
)
scK4me1w6.new.proj.rmdoub <- addHarmony(
  ArchRProj = scK4me1w6.new.proj.rmdoub,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  force = TRUE,
  groupBy = "Sample"
)
scK4me1w6.new.proj.rmdoub <- addClusters(
  input = scK4me1w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  force = TRUE,
  resolution = 1,
  maxClusters = 100,
  dimsToUse = 2:30
)
scK4me1w6.new.proj.rmdoub <- addUMAP(
  ArchRProj = scK4me1w6.new.proj.rmdoub,
  reducedDims = "Harmony",
  name = "IterativeLSIUMAP",
  force = TRUE,
  threads = 8,
  dimsToUse = 2:30,
  nNeighbors = 40,
  minDist = 0.7
)
p2ac <- plotEmbedding(
  ArchRProj = scK4me1w6.new.proj.rmdoub,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "IterativeLSIUMAP",
  log2Norm = FALSE,
  rastr = FALSE,
  size = 0.5,
  labelMeans = TRUE
)
p2ac
