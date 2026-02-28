
# Integrated R Script: NK Cell Subset Extraction and Data Format Conversion

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)


# Step 1: Extract and process NK cell subset
cat("Step 1: Extracting and processing NK cell subset...\n")

# Create output directory if it doesn't exist
dir.create("output/data", recursive = TRUE, showWarnings = FALSE)

# Read original data
Multi_myeloma <- readRDS("SigflowVelo/benchmark/GSE164551(MM)/data/TFVelo_data/GSE164551_seurat_afterAnno.rds")

# Extract NK cell subset
Multi_myeloma_NK <- subset(Multi_myeloma, subset = celltype %in% c("NK cells"))
cat(paste0("Extracted ", ncol(Multi_myeloma_NK), " NK cells\n"))

# Standard single-cell analysis pipeline
Multi_myeloma_NK <- FindVariableFeatures(
  Multi_myeloma_NK, 
  selection.method = "vst", 
  nfeatures = 3000, 
  assay = 'RNA'
)

Multi_myeloma_NK <- ScaleData(
  Multi_myeloma_NK, 
  verbose = FALSE, 
  assay = 'RNA'
)

Multi_myeloma_NK <- RunPCA(
  Multi_myeloma_NK, 
  npcs = 30, 
  verbose = FALSE
)

Multi_myeloma_NK <- RunUMAP(
  Multi_myeloma_NK, 
  reduction = "pca", 
  dims = 1:30,
  metric = "euclidean",
  n.epochs = 500,
  negative.sample.rate = 15L
)

Multi_myeloma_NK <- RunTSNE(
  Multi_myeloma_NK, 
  reduction = "pca", 
  dims = 1:30,
  negative.sample.rate = 15L
)

Multi_myeloma_NK <- FindNeighbors(
  Multi_myeloma_NK, 
  reduction = "pca", 
  dims = 1:30
)

Multi_myeloma_NK <- FindClusters(
  Multi_myeloma_NK, 
  resolution = 0.1, 
  verbose = FALSE
)

# Save processed NK cell object
saveRDS(Multi_myeloma_NK, file = "SigflowVelo/benchmark/GSE164551(MM)/data/TFVelo_data/Multi_myeloma_NK.rds")
cat("NK cell subset saved to: SigflowVelo/benchmark/GSE164551(MM)/data/TFVelo_data/Multi_myeloma_NK.rds\n")


# Step 2: Data format conversion (Seurat → h5ad)
cat("\nStep 2: Converting data format...\n")

# Set default assay to RNA
DefaultAssay(Multi_myeloma_NK) <- "RNA"

# Get raw count matrix
raw_counts <- GetAssayData(Multi_myeloma_NK, slot = "counts")

# Normalize data (if not already normalized)
Multi_myeloma_NK <- NormalizeData(Multi_myeloma_NK, verbose = FALSE)

# Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(Multi_myeloma_NK, assay = "RNA", data = "data")

# Add normalized counts and raw counts to assay
assays(sce)$counts <- logcounts(sce)  # Normalized counts
assays(sce)$count <- raw_counts       # Raw counts

# Save as h5ad format
output_path <- "SigflowVelo/benchmark/GSE164551(MM)/data/TFVelo_data/Multi_myeloma_NK.h5ad"
writeH5AD(
  sce,
  file = output_path,
  compression = "lzf"
)

cat(paste0("Data successfully converted to h5ad format and saved to: ", output_path, "\n"))
cat(paste0("Output file dimensions:\n"))
cat(paste0("  Rows (cells): ", ncol(sce), "\n"))
cat(paste0("  Columns (genes): ", nrow(sce), "\n"))