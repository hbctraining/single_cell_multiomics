# Demultiplexing HTO and protein expression ADT
# Single-cell multiomics modules
# Author: Harvard Chan Bioinformatics Core
# Created: October 2026

# Load libraries
library(tidyverse)
library(Seurat)

# Set seed
set.seed(123456789L)

# # Create directory tree
# fs::dir_tree("data/", recurse = 1)

# Load sample metadata
df_meta <- read_tsv("data/samples.tsv")
df_meta

# Load libraries and data
library(tidyverse)
library(Seurat)

# Set seed
set.seed(123456789L)
# Turn off scientific notation
options(scipen = 999)

# Load sample metadata
df_meta <- read_tsv("data/samples.tsv")

# Sample of interest
sample <- "PT1"

# Load raw counts data
counts <- Read10X(file.path("data", sample, "ADT"))
names(counts)

# RNA expression matrix
counts_gex <- counts[["Gene Expression"]]

# Create seurat object
seurat <- CreateSeuratObject(
  counts = counts_gex,
  assay = "RNA",
  project = sample)

seurat

# Protein expression matrix (ADT and HTO)
counts_adt <- counts[["Antibody Capture"]]

# View first few protein names
rownames(counts_adt) %>% head()

# Remove appended _TotalA string from protein names
rownames(counts_adt) <- str_replace_all(rownames(counts_adt), "_TotalA", "")

# Looking for hashes in the dataset
hashes <- rownames(counts_adt)[str_detect(rownames(counts_adt), "hashtag")]
hashes

# Map of hash to patient
df_hash <- read.csv(file.path("data",
                              sample,
                              "ADT", "hashing_info.txt.gz"),
                    header = FALSE, sep = ":",
                    col.names = c("hash", "patient")) %>%
  # data wrangling
  mutate(hash = str_replace(hash, "H", "h"),
         hash = str_replace(hash, " ", ""))
df_hash

# Grab counts for just the hashes used in this sample
counts_hash <- counts_adt[df_hash$hash, ]

# Add HTO assay to our seurat object
seurat[["HTO"]] <- CreateAssayObject(counts = counts_hash)
seurat

# Name of features in the HTO assay
Features(seurat[["HTO"]])

# Protein expression
counts_protein <- counts_adt[!rownames(counts_adt) %in% hashes, ]
seurat[["ADT"]] <- CreateAssayObject(counts = counts_protein)

seurat

# Set default assay to HTO
DefaultAssay(seurat) <- "HTO"

# Seurat metadata
seurat@meta.data %>% head()

ggplot(seurat@meta.data) +
  geom_histogram(aes(x = nCount_HTO),
                 fill = "darkred",
                 bins = 100) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() + NoLegend() +
  geom_vline(xintercept = 30)

# Filter based on HTO counts
seurat <- subset(seurat, nCount_HTO > 30)
# Filter based on RNA counts
seurat <- subset(seurat, nCount_RNA >= 500)
seurat <- subset(seurat, nCount_RNA <= 10000)
seurat

# CLR HTO normalization
seurat <- NormalizeData(
    seurat,
    assay = "HTO",
    normalization.method = "CLR",
    margin = 2
  )
seurat

# HTO demultiplexing
seurat <- HTODemux(
    seurat,
    assay = "HTO",
    positive.quantile = 0.99,
    verbose = TRUE)

# Seurat metadata (HTODemux columns only)
seurat@meta.data %>%
  select(starts_with("HTO_"), hash.ID) %>%
  head()

# Number of cells per demultiplexed hash ID
ggplot(seurat@meta.data) +
  geom_bar(aes(x = hash.ID, fill = hash.ID),
           color = "black") +
  # Add count above bars
  geom_text(aes(x = hash.ID, label = after_stat(count)),
          stat = "count", vjust = -0.5,
          size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Distribution of HTO counts per demultiplexed hash ID
VlnPlot(seurat,
        features = "nCount_HTO",
        group.by = "hash.ID",
        assay = "HTO",
        pt.size = 0) +
  scale_y_log10()

# HTO expression by demultiplexed hash ID
RidgePlot(
  seurat,
  features = df_hash$hash,
  group.by = "hash.ID",
  assay = "HTO")

# Remove doublets and negatives
seurat_demux <- subset(seurat, subset = (hash.ID %in% df_hash$hash))
seurat_demux

# Add patient to metadata by matching hash
meta <- left_join(seurat_demux@meta.data,
                  df_hash, by = c("hash.ID" = "hash"))
rownames(meta) <- Cells(seurat_demux)

seurat_demux@meta.data <- meta

# Load the processed seurat object
seurat_demux <- readRDS("data/processed_seurat_objects/seurat_demux.RDS")
seurat_demux

# Number of cells per patient
ggplot(seurat_demux@meta.data) +
  geom_bar(aes(x = patient,
               fill = hiv_status),
           color = "black") +
  geom_text(aes(x = patient,
                label = after_stat(count)),
          stat = "count", vjust = -0.5,
          size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Load libraries and data
library(Seurat)
library(tidyverse)
library(harmony)

# Set seed
set.seed(123456789L)

# Load demultiplexed Seurat object
seurat_demux <- readRDS("data/processed_seurat_objects/seurat_demux.RDS")

# Set default assay to ADT
DefaultAssay(seurat_demux) <- "ADT"

# Isotype control antibodies in the panel
isotypes <- c("ArmenianHamsterIgG",
              "Mouse-IgG1-isotype", "Mouse-IgG2a-isotype", "Mouse-IgG2b-isotype",
              "RatIgG1k", "RatIgG1l", "RatIgG2a", "Rat-IgG2b-isotype", "RatIgG2c")

# Isotype control expression by HIV status
VlnPlot(seurat_demux,
        group.by = "hiv_status",
        features = isotypes,
        assay = "ADT",
        pt.size = 0)

# CLR normalization across proteins within each cell
seurat_demux <- NormalizeData(
  seurat_demux,
  normalization.method = "CLR",
  margin = 2)

seurat_demux

# Use all proteins except the isotype controls
adt_features <- setdiff(rownames(seurat_demux[["ADT"]]), isotypes)
VariableFeatures(seurat_demux) <- adt_features

length(adt_features)

# Scale protein expression
seurat_demux <- ScaleData(seurat_demux)
seurat_demux

# PCA on protein expression
seurat_demux <- RunPCA(
  seurat_demux,
  npcs = 30,
  assay = "ADT",
  reduction.name = "pca_adt",
  reduction.key = "adtPC_",
  verbose = FALSE,
  seed.use = 123456789L)

# Protein PCA by library and HIV status
DimPlot(seurat_demux,
        group.by = c("orig.ident", "hiv_status"),
        reduction = "pca_adt",
        shuffle = TRUE)

# Top proteins for the first 5 PCs
print(seurat_demux[["pca_adt"]],
      dims = 1:5, nfeatures = 5)

# UMAP on uncorrected protein PCA
seurat_demux <- RunUMAP(
  seurat_demux,
  reduction = "pca_adt",
  dims = 1:20,
  reduction.name = "umap_adt_uncorrected",
  reduction.key = "adtUncorrectedUMAP_",
  verbose = FALSE,
  seed.use = 123456789L)

# Protein UMAP before integration by library and HIV status
DimPlot(seurat_demux,
        group.by = c("orig.ident", "hiv_status"),
        reduction = "umap_adt_uncorrected",
        shuffle = TRUE)

# Correct protein PCs for library
seurat_demux <- RunHarmony(
  seurat_demux,
  group.by.vars = "orig.ident",
  reduction.use = "pca_adt",
  dims.use = 1:20,
  reduction.save = "harmony_adt",
  verbose = FALSE)

# UMAP on corrected protein PCs
seurat_demux <- RunUMAP(
  seurat_demux,
  reduction = "harmony_adt",
  dims = 1:20,
  reduction.name = "umap_adt",
  reduction.key = "adtUMAP_",
  verbose = FALSE,
  seed.use = 123456789L)

# Protein UMAP after integration by library and HIV status
DimPlot(seurat_demux,
        group.by = c("orig.ident", "hiv_status"),
        reduction = "umap_adt",
        shuffle = TRUE)

# Build the nearest neighbor graph on corrected protein PCs
seurat_demux <- FindNeighbors(
  seurat_demux,
  reduction = "harmony_adt",
  dims = 1:20,
  graph.name = c("nn_adt", "snn_adt"),
  verbose = FALSE)

# Cluster the cells at a resolution of 0.8
seurat_demux <- FindClusters(
  seurat_demux,
  graph.name = "snn_adt",
  resolution = 0.8,
  random.seed = 123456789L,
  verbose = FALSE)

# RNA clusters on the RNA UMAP
p_rna <- DimPlot(seurat_demux,
                 group.by = "snn_rna_res.0.8",
                 reduction = "umap_rna",
                 label = TRUE) +
  ggtitle("RNA clusters")

# Protein clusters on the protein UMAP
p_adt <- DimPlot(seurat_demux,
                 group.by = "snn_adt_res.0.8",
                 reduction = "umap_adt",
                 label = TRUE) +
  ggtitle("Protein clusters")

# Plot side by side
p_rna + p_adt

# RNA clusters on the protein UMAP
DimPlot(seurat_demux,
        group.by = "snn_rna_res.0.8",
        reduction = "umap_adt",
        label = TRUE)

# CD4 and CD8 expression from RNA and protein
FeaturePlot(seurat_demux,
            features = c("rna_CD4", "rna_CD8A", "adt_CD4", "adt_CD8"),
            reduction = "umap_rna",
            order = TRUE,
            ncol = 2)

# Marker genes for each cell type
markers_gut <- list(
  "CD4 T"                 = c("CD3E", "CD4", "IL7R", "CD40LG"),
  "CD8 T"                 = c("CD8A", "CD8B", "GZMK"),
  "B naive"               = c("MS4A1", "IGHD", "TCL1A", "FCER2"),
  "B memory"              = c("CD27", "TNFRSF13B", "AIM2"),
  "B proliferating"       = c("MKI67", "TOP2A", "AICDA", "RGS13"),
  "Plasma cell - IgA"     = c("JCHAIN", "IGHA1", "IGHA2", "MZB1"),
  "Plasma cell - IgG"     = c("IGHG1", "IGHG3", "IGHG4"),
  "CD14 monocytes"        = c("CD14", "LYZ", "S100A8", "FCN1"),
  "NK"                    = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
  "Mast"                  = c("TPSAB1", "CPA3", "KIT", "MS4A2"),
  "Epithelial"            = c("EPCAM", "KRT8", "KRT20", "CDH1"),
  "Colonic subepithelial" = c("PDGFRA", "F3", "BMP4", "SOX6"),
  "Endothelial"           = c("PECAM1", "VWF", "CDH5", "PLVAP"),
  "Enteric glia"          = c("PLP1", "S100B", "SOX10", "CDH19"),
  "Myenteric ganglia"     = c("ELAVL4", "PHOX2B", "RET", "SNAP25"),
  "Myofibroblast"         = c("ACTA2", "TAGLN", "MYH11", "DES"))

DotPlot(seurat_demux, markers_gut,
        group.by = "snn_rna_res.0.8",
        assay = "RNA") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Cluster -> cell type labels
celltype_labels <- list(
  "1"  = "CD4 T",
  "2"  = "Plasma cell - IgA",
  "3"  = "CD8 T",
  "4"  = "B memory",
  "5"  = "Plasma cell - IgA",
  "6"  = "NK",
  "7"  = "Plasma cell - IgA",
  "8"  = "Plasma cell - IgA",
  "9"  = "Colonic subepithelial",
  "10" = "B naive",
  "11" = "CD4 T",
  "12" = "Plasma cell - IgG",
  "13" = "B proliferating",
  "14" = "Colonic subepithelial",
  "15" = "CD14 monocytes",
  "16" = "Plasma cell - IgA",
  "17" = "Endothelial",
  "18" = "Enteric glia",
  "19" = "Mast",
  "20" = "Myofibroblast"
)

# Color for each cell type
pal_gut <- c(
  "CD4 T"                 = "#F4766E",
  "CD8 T"                 = "#E78526",
  "B naive"               = "#FCDE13",
  "B memory"              = "#D8E025",
  "B proliferating"       = "#8FC740",
  "Plasma cell - IgA"     = "#0C964A",
  "Plasma cell - IgG"     = "#106B39",
  "CD14 monocytes"        = "#1CBAE9",
  "NK"                    = "#1374BC",
  "Mast"                  = "#273B90",
  "Epithelial"            = "#808EC5",
  "Colonic subepithelial" = "#A780BB",
  "Endothelial"           = "#E16CA8",
  "Enteric glia"          = "#C79A6B",
  "Myenteric ganglia"     = "#AB7D52",
  "Myofibroblast"         = "#754D29"
)

# Set idents to the RNA clusters
Idents(seurat_demux) <- "snn_rna_res.0.8"

# Rename the clusters and keep the labels in the metadata
seurat_demux <- RenameIdents(seurat_demux, celltype_labels)
seurat_demux$celltype <- Idents(seurat_demux)
seurat_demux$celltype <- factor(seurat_demux$celltype, levels = names(pal_gut))

# Count number of cells
table(seurat_demux$celltype)

# Cell types on the RNA UMAP
p1 <- DimPlot(seurat_demux,
        group.by = "celltype",
        reduction = "umap_rna",
        cols = pal_gut) + 
  NoLegend() +
  ggtitle("UMAP RNA")

# Cell types on the protein UMAP
p2 <- DimPlot(seurat_demux,
        group.by = "celltype",
        reduction = "umap_adt",
        cols = pal_gut) +
  ggtitle("UMAP ADT")

p1 | p2

# RNA and protein expression of CD4 for each cell
df_cd4 <- FetchData(seurat_demux,
                    vars = c("rna_CD4", "adt_CD4", "celltype"))
head(df_cd4)

df_cd4 <- subset(df_cd4, celltype == "CD4 T")

# Compare CD4 RNA and protein for CD4+ T cells
ggplot(df_cd4) +
  geom_point(aes(x = rna_CD4,
                 y = adt_CD4,
                 color = celltype),
             size = 0.5,
             alpha = 0.5) +
  scale_color_manual(values = pal_gut) +
  labs(x = "CD4 RNA expression",
       y = "CD4 protein expression") +
  theme_bw() +
  NoLegend()
