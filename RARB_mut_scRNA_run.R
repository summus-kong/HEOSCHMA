library(Seurat)
library(harmony)

dir <- "./"
dirs <- dir(dir, "[1-9]$")
sample_dirs <- paste(paste0(dir, dirs), "outs/filtered_feature_bc_matrix", sep = "/")

sclist <- lapply(dirs, function(x) {
  print(x)
  sce <- CreateSeuratObject(
    counts = Read10X(paste(
      paste0(dir, x), "outs/filtered_feature_bc_matrix", sep = "/"
    )),
    project = str_replace(x, "_", ""),
    min.cells = 3,
    min.features = 200
  )
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce[["percent.hb"]] <- PercentageFeatureSet(sce, pattern = "^HB[^(P)]")
  return(sce)
})
names(sclist) <- str_replace(dirs, "_", "")
sclist[[1]]
rownames(sclist[[1]])
summary(sclist[[1]]@meta.data$nFeature_RNA)

rarb_scrna <- merge(sclist[[1]], sclist[2], add.cell.ids = str_replace(dirs, "_", ""))

feats <- c("nCount_RNA",
           "nFeature_RNA",
           "percent.mt",
           "percent.ribo",
           "percent.hb")
VlnPlot(
  rarb_scrna,
  features = feats,
  pt.size = 0,
  group.by = "orig.ident",
  ncol = 3
)
rarb_scrna1 <- subset(
  rarb_scrna,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 8000 &
    percent.mt < 30 &
    percent.ribo < 25 &
    nCount_RNA < quantile(nCount_RNA, 0.97) &
    nCount_RNA > 600
)

table(rarb_scrna1$orig.ident)
VlnPlot(
  rarb_scrna1,
  features = feats,
  pt.size = 0,
  group.by = "orig.ident",
  ncol = 3
)

library(DoubletFinder)
scRNAlist <- SplitObject(rarb_scrna1, split.by = "orig.ident")
scRNAlist <- lapply(
  X = scRNAlist,
  FUN = function(pbmc) {
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc)
    pbmc <- RunUMAP(pbmc, dims = 1:30)
    pbmc <- FindNeighbors(pbmc, dims = 1:30)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    sweep.res.list_pbmc <- paramSweep_v3(pbmc, PCs = 1:30, sct = FALSE)
    sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
    mpK <- as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
    annotations <- pbmc$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(0.015 * nrow(pbmc@meta.data))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    pbmc <- doubletFinder_v3(
      pbmc,
      PCs = 1:30,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi,
      reuse.pANN = FALSE
    )
    pbmc@meta.data[, "DF_hi.lo"] <- pbmc@meta.data[, ncol(pbmc@meta.data)]
    return(pbmc)
  }
)

for (i in 1:length(scRNAlist)) {
  print(i)
  scRNAlist[[i]] <- SCTransform(scRNAlist[[i]])
}
features1 <- SelectIntegrationFeatures(object.list = scRNAlist)
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = features1)

sclist.anchors <- FindIntegrationAnchors(
  object.list = scRNAlist,
  dims = 1:30,
  normalization.method = "SCT",
  reduction = "cca",
  anchor.features = features1
)
rarb65893_scrna.cca <- IntegrateData(
  anchorset = sclist.anchors,
  new.assay.name = "CCA",
  normalization.method = "SCT",
  dims = 1:30
)
rarb65893_scrna.cca1 <- rarb65893_scrna.cca %>%
  RunPCA()
rarb65893_scrna.cca1 <- rarb65893_scrna.cca1 %>%
  RunUMAP(dims = 1:20) %>%
  RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)

DefaultAssay(rarb65893_scrna.cca1) <- "RNA"
rarb65893_scrna.cca1 <- NormalizeData(rarb65893_scrna.cca1) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(.))

DimPlot(
  rarb65893_scrna.cca1,
  reduction = "umap",
  label = TRUE,
  group.by = "seurat_clusters"
) +
  DimPlot(rarb65893_scrna.cca1,
          reduction = "umap",
          group.by = "orig.ident") +
  DimPlot(rarb65893_scrna.cca1,
          reduction = "umap",
          group.by = "DF_hi.lo")

identitymappings <- c(
  "0" = "EC",
  "1" = "fibroblast",
  "2" = "CMs",
  "3" = "fibroblast",
  "4" = "fibroblast",
  "5" = "Peric",
  "6" = "CMs",
  "7" = "EC",
  "8" = "EndoC",
  "10" = "macrophages",
  "11" = "blood"
)
rarb65893_scrna.cca1 <- RenameIdents(rarb65893_scrna.cca1, identitymappings)
rarb65893_scrna.cca1$CellType <- Idents(rarb65893_scrna.cca1)
DefaultAssay(rarb65893_scrna.cca1) <- "CCA"
rarb65893_scrna.cca1 <- rarb65893_scrna.cca1 %>%
  RunPCA()
rarb65893_scrna.cca1 <- rarb65893_scrna.cca1 %>%
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5)

DefaultAssay(rarb65893_scrna.cca1) <- "RNA"
rarb65893_scrna.cca1 <- NormalizeData(rarb65893_scrna.cca1) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(.))
colsab <- c("#DC5153",
            "#396AA8",
            "#EED01E",
            "#62A9A1",
            "#4AA259",
            "#F2C1BC",
            "#DE6830")
names(colsab) <- c("CMs",
                   "EC",
                   "fibroblast",
                   "Peric",
                   "EndoC",
                   "macrophages",
                   "blood")
unique(rarb65893_scrna.cca1$CellType)
DimPlot(
  rarb65893_scrna.cca1,
  reduction = "umap",
  label = FALSE,
  group.by = "CellType",
  cols = colsab
)
DimPlot(
  rarb65893_scrna.cca1,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.3
)
