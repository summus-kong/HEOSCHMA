options(scipen = 100)

library(ArchR)
library(ArchRtoSignac)
library(Pando)
library(Seurat)
library(tidyverse)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
library(EnsDb.Hsapiens.v86)

registerDoParallel(4)

data("human_pwms_v2")
data("phastConsElements20Mammals.UCSC.hg38")

ArchR2Signac_us <- function(ArchRProject,
                            refversion,
                            samples = NULL,
                            fragments_dir = NULL,
                            pm,
                            fragments_fromcellranger = NULL,
                            fragments_file_extension = NULL,
                            keyvalue = NULL,
                            annotation) {
  if (is.null(samples)) {
    samples <- unique(ArchRProject@cellColData$Sample)
  }

  if (fragments_fromcellranger == "YES" |
      fragments_fromcellranger == "Y" |
      fragments_fromcellranger == "Yes") {
    print("In Progress:")
    print("Prepare Seurat list for each sample")

    output_dir = "/outs/"

    seurat_list <- lapply(samples, function(cur_sample) {
      print(cur_sample)
      cur_fragments <- ifelse(
        is.list(fragments_dir),
        paste0(fragments_dir[[which(samples == cur_sample)]], output_dir, "fragments.tsv.gz"),
        paste0(
          fragments_dir,
          cur_sample,
          output_dir,
          "fragments.tsv.gz"
        )
      )

      cur_pm <- pm[, grepl(paste0("^", cur_sample, "#"), colnames(pm))]
      cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)

      colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), "#"))[, 2]
      rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), "#"))[, 2]
      print(dim(cur_pm))
      cur_chromatin <- Signac::CreateChromatinAssay(
        counts = cur_pm,
        sep = c("-", "-"),
        fragments = cur_fragments,
        ranges = ArchRProject@peakSet,
        genome = refversion,
        annotation = annotation
      )

      cur_atac <- Seurat::CreateSeuratObject(cur_chromatin, assay = "peaks", meta.data = cur_meta, )

    })
  }

  if (fragments_fromcellranger == "NO" |
      fragments_fromcellranger == "N" |
      fragments_fromcellranger == "No") {
    print("IF selecting NO, please make sure to provide fragments_file_extension")
    print("In Progress:")
    print("Prepare Seurat list for each sample")

    seurat_list <- lapply(samples, function(cur_sample) {
      print(cur_sample)
      cur_sample_value <- basename(keyvalue[cur_sample])
      print(cur_sample_value)
      cur_fragments <- ifelse(
        is.list(fragments_dir),
        paste0(fragments_dir[[which(samples == cur_sample)]], fragments_file_extension),
        paste0(
          fragments_dir,
          cur_sample_value,
          fragments_file_extension
        )
      )
      print(cur_fragments)
      cur_pm <- pm[, grepl(paste0("^", cur_sample, "#"), colnames(pm))]
      cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)

      colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), "#"))[, 2]
      rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), "#"))[, 2]
      print(dim(cur_pm))
      cur_chromatin <- Signac::CreateChromatinAssay(
        counts = cur_pm,
        sep = c("-", "-"),
        ranges = ArchRProject@peakSet,
        genome = refversion,
        annotation = annotation
      )

      cur_atac <- Seurat::CreateSeuratObject(cur_chromatin, assay = "peaks", meta.data = cur_meta, )

    })

  }


  print("In Progress:")
  print("Merge Seurat list")
  SeuratObject <- if (length(seurat_list) > 1) {
    merge(x = seurat_list[[1]],
          y = seurat_list[2:length(seurat_list)],
          add.cell.ids = samples)
  } else {
    seurat_list[[1]]
  }

  print("Return SeuratObject")
  SeuratObject

}

annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

scK27w46.1.proj_merge.rmdoub.card
scK27w46.1.proj_merge.rmdoub.card <- addImputeWeights(scK27w46.1.proj_merge.rmdoub.card)
pkm <- getPeakMatrix(scK27w46.1.proj_merge.rmdoub.card)

fragments_dir <- "./"
inputFiles <- c(
  "../sck27ac-w6-e1-fragmentsInCells.tsv.gz",
  "../sck27ac-w6-e2-fragmentsInCells.tsv.gz",
  "../sck27ac-w6-e3-fragmentsInCells.tsv.gz",
  "sck27ac-w4-e1-fragmentsInCells.tsv.gz",
  "sck27ac-w4-e2-fragmentsInCells.tsv.gz",
  "sck27ac-w4-e3-fragmentsInCells.tsv.gz"
)
names(inputFiles) <- c(
  "humanw6-1me",
  "humanw6-2me",
  "humanw6-3me",
  "humanw4-k27ac-1me",
  "humanw4-k27ac-2me",
  "humanw4-k27ac-3me"
)
seurat_atac_card2507 <- ArchR2Signac_us(
  ArchRProject = scK27w46.1.proj_merge.rmdoub.card,
  refversion = "hg38",
  fragments_dir = fragments_dir,
  pm = pkm,
  fragments_fromcellranger = "NO",
  fragments_file_extension = NULL,
  keyvalue = inputFiles,
  annotation = annotations
)

gsm.card2507 <- getGeneScoreMatrix_us(ArchRProject = scK27w46.1.proj_merge.rmdoub.card,
                                      SeuratObject = seurat_atac_card2507,
                                      useMatrix = "GeneIntegrationMatrix")
dim(gsm.card2507)
seurat_atac_card2507[["RNA"]] <- CreateAssayObject(counts = gsm.card2507)

umap_df <- scK27w46.1.proj_merge.rmdoub.card@embeddings[["IterativeLSIUMAP"]]$df %>% as.matrix
rownames(umap_df) <- colnames(seurat_atac_card2507)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
seurat_atac_card2507@reductions$umap <- Seurat::CreateDimReducObject(embeddings = umap_df, assay = "peaks")
LSI_matrix  <- scK27w46.1.proj_merge.rmdoub.card@reducedDims[["IterativeLSI"]][[1]]
LSI_matrix <- LSI_matrix[rownames(LSI_matrix) %in% str_replace_all(colnames(seurat_atac_card2507), "_", "#"), ]
rownames(LSI_matrix) <- colnames(seurat_atac_card2507)
colnames(LSI_matrix) <- paste0("LSI_", 1:ncol(LSI_matrix))
seurat_atac_card2507@reductions$IterativeLSI <- Seurat::CreateDimReducObject(embeddings = LSI_matrix, assay = "peaks")
data("phastConsElements20Mammals.UCSC.hg38")
seurat_atac_card2507 <- initiate_grn(
  seurat_atac_card2507,
  rna_assay = "RNA",
  peak_assay = "peaks",
  regions = phastConsElements20Mammals.UCSC.hg38
)
seurat_atac_card2507
GetGRN(seurat_atac_card2507)
regions <- NetworkRegions(seurat_atac_card2507)
regions@ranges

seurat_atac_card2507 <- find_motifs(seurat_atac_card2507, pfm = human_pwms_v2, genome = BSgenome.Hsapiens.UCSC.hg38)



seurat_atac_card2507 <- infer_grn(
  seurat_atac_card2507,
  peak_to_gene_method = "GREAT",
  downstream = 1e+05,
  upstream = 1e+05,
  tf_cor = 0.03,
  method = "glmnet",
  parallel = TRUE
)
coef(seurat_atac_card2507)

seurat_atac_card2507 <- find_modules(
  seurat_atac_card2507,
  p_thresh = 0.05,
  min_genes_per_module = 1,
  rsq_thresh = 0.1
)

seurat_atac_card25071 <- get_network_graph(
  seurat_atac_card2507,
  rna_assay = "RNA",
  n_neighbors = 2,
  features = cm_features,
  umap_method = "corr"
)
subgraph.card <- NetworkGraph(seurat_atac_card25071)
subgraph.card
igraph::as_long_data_frame(subgraph.card) -> v.card
library(ggraph)
ggraph(subgraph.card, layout = "drl") +
  geom_edge_parallel(aes(color = corr, alpha = after_stat(index)),
                     width = 0.8,
                     sep = unit(0.8, "mm")) +
  geom_node_point(
    aes(size = centrality),
    color = "black",
    fill = "#778899",
    shape = 21,
    stroke = 0.2
  ) +
  geom_node_text(aes(label = name),
                 size = 4,
                 shape = 21,
                 repel = TRUE) +
  scale_edge_colour_gradient2(low = "grey",
                              high = "#cac0e1",
                              mid = "#f3dac0") +
  theme_test()
