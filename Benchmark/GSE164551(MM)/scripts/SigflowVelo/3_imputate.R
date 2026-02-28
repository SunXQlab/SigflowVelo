# Load required packages
library(Seurat)
library(scImpute)

# Read count data from previous step
scimpute(
  # Full path to raw count matrix
  count_path = 'SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data/count_file.csv', 
  infile = "csv",           # Format of input file
  outfile = "csv",          # Format of output file
  out_dir = "SigflowVelo/benchmark/GSE164551(MM)/data/SigflowVelo_data", # Full path to output directory
  labeled = FALSE,          # Cell type labels not available
  drop_thre = 0.5,          # Threshold set on dropout probability
  Kcluster = 2,             # Number of cell subpopulations
  ncores = 1)               # Number of cores to use

