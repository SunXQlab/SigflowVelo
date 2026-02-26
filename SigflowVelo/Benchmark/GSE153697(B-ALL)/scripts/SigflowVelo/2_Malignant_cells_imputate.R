library(Seurat)
library(scImpute)

scimpute(# full path to raw count matrix
  count_path = "../../resluts/SigflowVelo/adata/intermediate_results/count_file.csv", 
  infile = "csv",           # format of input file
  outfile = "csv",          # format of output file
  out_dir = "../../resluts/SigflowVelo/adata/intermediate_results/scImputeOutput/",     #path to output directory
  labeled = FALSE,          # cell type labels not available
  drop_thre = 0.5,          # threshold set on dropout probability
  Kcluster = 2,             # 2 cell subpopulations
  ncores = 1) 
