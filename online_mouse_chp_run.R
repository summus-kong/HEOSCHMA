

library(tidyverse)
library(harmony)
library(Seurat)



chp_gse168704 <- CreateSeuratObject(
  counts = Read10X("./chp/gse168704/scrna/", gene.column = 1),
  min.cells = 3,
  min.features = 200,
  project = "chp"
)
chp_meta <- read_table(gzfile("./chp/gse168704/scrna/GSE168704_sc_metadata.tsv.gz"),
                       col_names = TRUE) %>%
  column_to_rownames(var = "cellBC") %>%
  dplyr::slice(match(rownames(chp_gse168704@meta.data), rownames(.))) %>%
  as.data.frame()


chp_gse168704@meta.data <- cbind(chp_gse168704@meta.data, chp_meta)
chp_gse168704 %<>%
  NormalizeData() %>%
  FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  harmony::RunHarmony(., group.by.vars = "orig.ident") %>%
  RunUMAP(., dims = 1:20, reduction = "harmony") %>%
  FindNeighbors(., dims = 1:20, reduction = "harmony") %>%
  FindClusters(., resolution = 0.2)

chp_pals <- c(
  "Epithelial" = "#183D68",
  "Mesenchymal" = "#b2db87",
  "NeuroGlia" = "#e187cb",
  "Endothelial" = "#b24646",
  "Immune" = "#64cccf"
)

DimPlot(chp_gse168704, group.by = "cell.id", cols = chp_pals) +
  labs(title = "GSE168704")

tfs <- c("OTX2", "TEAD1")
FeaturePlot(
  chp_gse168704,
  order = FALSE,
  features = str_to_title(tfs),
  reduction = "umap",
  ncol = 3
)

library(Nebulosa)
library(ggrastr)
Nebulosa::plot_density(
  chp_gse168704,
  features = "Otx2",
  reduction = "umap",
  pal = "cividis",
  size = 0.5,
  na.value = "#F6EEB3"
)

plist <- list()
for (i in seq_along(str_to_title(tfs))) {
  p <- Nebulosa::plot_density(
    chp_gse168704,
    features = str_to_title(tfs)[i],
    pal = "cividis",
    raster = TRUE,
    size = 0.8
  ) &
    theme(legend.position = "none")

  plist[[i]] <- p
}
library(cowplot)
plot_grid(plotlist = plist, ncol = 4)

Idents(chp_gse168704) <- chp_gse168704$cell.id
chp_gse168704.markers <- FindAllMarkers(
  chp_gse168704,
  only.pos = TRUE,
  min.pct = 0.1,
  return.thresh = 0.05,
  logfc.threshold = 0.25
)

table(chp_gse168704$cell.id) %>%
  as.data.frame() %>%
  mutate(frac = Freq / sum(Freq)) %>%
  mutate(Var1 = fct_reorder(Var1, frac, .desc = TRUE)) %>%
  ggplot(aes(Var1, frac, fill = Var1)) +
  geom_col(color = "black", width = 0.7) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  ggsci::scale_fill_jama() +
  labs(x = NULL, y = "Cell Fraction") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.line.x = element_line(color = "black", size = 0.6),
    axis.line.y = element_line(color = "black", size = 0.6),
    axis.title = element_text(color = "black", size = 12),
    axis.text.x = element_text(
      color = "black",
      size = 10,
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "none"
  )
