options(scipen = 100)

library(Seurat)
library(Signac)
library(ArchR)
library(tidyverse)

pcw11_dir <- "/public/ojsys/brain/lingshukuan/houwb/cell_scCUT_revise/pcw11_mailuocong_scrnaatac"

pcw11_rna_atac <- readRDS(paste0(pcw11_dir, "/W11.rds"))
DefaultAssay(pcw11_rna_atac) <- "RNA"
pcw11_rna_atac[["percent.mt"]] <- PercentageFeatureSet(pcw11_rna_atac, pattern = "^MT-")
pcw11_rna_atac[["percent.ribo"]] <- PercentageFeatureSet(pcw11_rna_atac, pattern = "^RP[SL]")
VlnPlot(
  pcw11_rna_atac,
  group.by = "orig.ident",
  ncol = 3,
  features = c(
    "percent.mt",
    "nCount_RNA",
    "nFeature_RNA",
    "nCount_ATAC",
    "nFeature_ATAC"
  ),
  pt.size = 0,
  cols = "red"
)

DimPlot(pcw11_rna_atac, reduction = "umap")
FeaturePlot(
  pcw11_rna_atac,
  reduction = "umap",
  features = c(
    "HTR2C",
    "ENPP2",
    "KCNJ13",
    "AQP1",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "ADAM12",
    "DCN",
    "LUM",
    "PEAR1",
    "CXCL12",
    "FLT1",
    "PLVAP",
    "PECAM1",
    "PTPRC",
    "HMHA1",
    "P2RY6",
    "LY86",
    "SLC1A3"
  ),
  ncol = 5
)
pcw11_rna_atac1 <- subset(
  pcw11_rna_atac,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 10000 &
    percent.mt < 40 &
    percent.ribo < 25 &
    nCount_RNA < 40000 & nCount_RNA > 200
)

pcw11_rna_atac1$nCount_RNA %>%
  as.data.frame() %>%
  rename("nCount_RNA" = ".") %>%
  ggplot(aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 250,
                 fill = "#00DFFF",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(nCount_RNA)),
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  labs(title = "Number of UMI count per cell", x = "UMI count", y = "Cell frequency") +
  annotate(
    "text",
    x = mean(pcw11_rna_atac1$nCount_RNA),
    y = 500,
    label = paste("mean = ", round(mean(
      pcw11_rna_atac1$nCount_RNA
    ))),
    color = "red",
    hjust = -0.1,
    size = 5
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line.x = element_line(color = "black", size = 0.6),
    axis.line.y = element_line(color = "black", size = 0.6),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10)
  )

pcw11_rna_atac$nFeature_RNA %>%
  as.data.frame() %>%
  rename("nFeature_RNA" = ".") %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 150,
                 fill = "#00DFFF",
                 color = "black") +
  geom_vline(
    aes(xintercept = mean(nFeature_RNA)),
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  labs(title = "Number of UMI count per cell", x = "Gene count", y = "Cell frequency") +
  annotate(
    "text",
    x = mean(pcw11_rna_atac1$nFeature_RNA),
    y = 450,
    label = paste("mean = ", round(mean(
      pcw11_rna_atac1$nFeature_RNA
    ))),
    color = "red",
    hjust = -0.1,
    size = 5
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line.x = element_line(color = "black", size = 0.6),
    axis.line.y = element_line(color = "black", size = 0.6),
    axis.title = element_text(color = "black", size = 12),
    axis.text = element_text(color = "black", size = 10)
  )

identitymappings <- c(
  "0" = "Epithelial",
  "1" = "Mesenchymal",
  "2" = "Epithelial",
  "3" = "Epithelial",
  "4" = "Endothelial",
  "5" = "Immune",
  "6" = "Epithelial",
  "7" = "Neuronal/Glial",
  "8" = "Mesenchymal"
)
pcw11_rna_atac <- RenameIdents(pcw11_rna_atac, identitymappings)
pcw11_rna_atac$CellType <- Idents(pcw11_rna_atac)

chp_pals1 <- c(
  "Epithelial" = "#183D68",
  "Mesenchymal" = "#b2db87",
  "Neuronal/Glial" = "#e187cb",
  "Endothelial" = "#b24646",
  "Immune" = "#64cccf"
)

DimPlot(
  pcw11_rna_atac1,
  reduction = "umap",
  label = FALSE,
  group.by = "CellType",
  cols = chp_pals1
) +
  DimPlot(
    pcw11_rna_atac1,
    reduction = "atacumap",
    label = FALSE,
    group.by = "CellType",
    cols = chp_pals1,
    order = TRUE
  )

Idents(pcw11_rna_atac) <- pcw11_rna_atac$CellType
DefaultAssay(pcw11_rna_atac) <- "RNA"
pcw11_rna_atac.markers <- FindAllMarkers(
  pcw11_rna_atac,
  only.pos = TRUE,
  min.pct = 0.2,
  return.thresh = 0.05,
  logfc.threshold = 0.25
)

top50 <- pcw11_rna_atac.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup()
avg_exp <- AverageExpression(
  object = pcw11_rna_atac,
  features = unique(top50$gene),
  group.by = "CellType",
  assays = "RNA",
  slot = "data"
)

avg_matrix <- avg_exp$RNA
my_colors <- colorRampPalette(c("white", "#FFE5B4", "orange", "darkred"))(100)
library(pheatmap)
pheatmap(
  avg_matrix,
  scale = "row",
  color = my_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 8,
  main = "Top 50 Markers Average Expression",
  border_color = NULL,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

library(Nebulosa)
library(ggrastr)
Nebulosa::plot_density(
  pcw11_rna_atac,
  features = "OTX2",
  reduction = "umap",
  pal = "inferno",
  size = 0.5,
  na.value = "#F6EEB3"
)
tfs <- c("OTX2", "TEAD1")
plist <- list()
for (i in seq_along(tfs)) {
  p <- Nebulosa::plot_density(
    pcw11_rna_atac,
    features = tfs[i],
    reduction = "wnn.umap",
    pal = "inferno",
    raster = TRUE,
    size = 0.8,
    na.value = "#F6EEB3"
  ) &
    theme(legend.position = "none")

  plist[[i]] <- p
}
library(cowplot)
plot_grid(plotlist = plist, ncol = 4)

pcw11_rna_atac <- FindMultiModalNeighbors(
  pcw11_rna_atac,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50)
)
pcw11_rna_atac <- RunUMAP(
  pcw11_rna_atac,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)
DimPlot(
  pcw11_rna_atac,
  reduction = "umap",
  label = FALSE,
  group.by = "CellType",
  cols = chp_pals1
) +
  DimPlot(
    pcw11_rna_atac,
    reduction = "atacumap",
    label = FALSE,
    group.by = "CellType",
    cols = chp_pals1
  ) +
  DimPlot(
    pcw11_rna_atac,
    reduction = "wnn.umap",
    group.by = "CellType",
    label = FALSE,
    cols = chp_pals1,
    repel = TRUE
  ) + ggtitle("WNN")
