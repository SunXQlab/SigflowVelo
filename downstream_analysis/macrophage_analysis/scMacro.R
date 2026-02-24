rm(list = ls())
library(Seurat)
library(paletteer)
library(ggplot2)
library(dplyr)

library(tidydr)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
## Fig2

scPDAC <- readRDS("./Data/scPDAC_anno.rds")
# 
scPDAC_macro <- subset(scPDAC, subset = celltype %in% c("Macrophage"))
scPDAC_macro <- FindVariableFeatures(scPDAC_macro, selection.method = "vst", nfeatures = 3000,assay = 'RNA')
scPDAC_macro <- ScaleData(scPDAC_macro, verbose = FALSE,Dassay = 'RNA')
scPDAC_macro <- RunPCA(scPDAC_macro, npcs = 30, verbose = FALSE)
scPDAC_macro <- RunUMAP(scPDAC_macro,
                        reduction = "pca",
                        dims = 1:30,
                        metric = "euclidean", #manhattan
                        n.epochs = 500,
                        negative.sample.rate = 15L)
scPDAC_macro <- FindNeighbors(scPDAC_macro, reduction = "pca", dims = 1:30)
scPDAC_macro <- FindClusters(scPDAC_macro, resolution = 0.2,verbose = FALSE)
p.1 <- DimPlot(scPDAC_macro, reduction = "umap", label = TRUE)
ggsave("./Result/Macrophage/umap_myeloid.pdf", plot = p.1,width = 8,height = 8)

# Myeloid注释
macro_annotate_dotplot <- function(SeuratObject){
  features <- c("THBS1", "EREG", "S100A9", "FCN1", # Monocyte (FCN1)
                "AREG", "HLA-DPB1", "HLA-DQA2", "FCER1A", "CD1C",  "CST7", # DC (CD1C)
                "CCL4", "CCL4L2", "CCL3", "CCL3L1", # M1-like(IL6)
                "RNASE1", "F13A1", "CD68", "TGFB1", "CD163", "ARG1" # M2-like (MS4A6A)
  )
  
  # DotPlot
  DefaultAssay(SeuratObject)<-"RNA"
  p <- DotPlot(SeuratObject, features = features, group.by = "seurat_clusters") + 
    theme_bw(base_size = 13) + 
    theme(panel.grid = element_blank(),  
          axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) + #
    scale_color_gradientn(values = seq(0,1,0.2),
                          colours = c(alpha("#6593CE", 0.2), alpha("#70C9EB", 0.2),"#EDA7A7","#FF9896FF")) + 
    labs(x=NULL,y=NULL) 
  ggsave("./Result/Macrophage/anno_dotplot_2.pdf", plot = p,width = 12,height = 6)
}

macro_annotate_dotplot(scPDAC_macro)

DefaultAssay(scPDAC_macro)<-"RNA"
anno <- data.frame(  
  seurat_clusters = scPDAC_macro@meta.data$seurat_clusters,  
  celltype = NA 
)  
anno$celltype[anno$seurat_clusters %in% c(0)] <- "Macrophage"  
anno$celltype[anno$seurat_clusters %in% c(1)] <- "Macrophage"  
anno$celltype[anno$seurat_clusters %in% c(2)] <- "Macrophage"  
anno$celltype[anno$seurat_clusters %in% c(3)] <- "Monocyte"  
anno$celltype[anno$seurat_clusters %in% c(4)] <- "Macrophage"  
anno$celltype[anno$seurat_clusters %in% c(5)] <- "DC" 
anno$celltype[anno$seurat_clusters %in% c(6)] <- "undefined" 
anno$celltype[anno$seurat_clusters %in% c(7)] <- "undefined" 
anno$celltype[anno$seurat_clusters %in% c(8)] <- "undefined" 

scPDAC_macro@meta.data$celltype2 <- anno$celltype

#========== umap myeloid ==============#

pp6 <- DimPlot(scPDAC_macro, label = F, group.by = "celltype2")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Macrophage/umap_myeloid.pdf',plot = pp6, width = 4, height = 4)


scMacro <- readRDS("./Data/Macrophage/scMacro.rds")

#========== umap group ==============#
my_col <- c("#E95D69", "#5CACEE")

pp4 <- DimPlot(scMacro, label = F, group.by = "group")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Macrophage/umap_group.pdf',plot = pp4, width = 4, height = 4)

#============ umap anno ==============#
my_col <- c("#E26CA8", "#8D78B7", "#E5007E", "#6593CE", "#31B6A9", "#70C9EB")
pp5 <- DimPlot(scMacro, label = T, group.by = "anno", label.size = 5) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=my_col)
ggsave('./Result/Macrophage/umap_anno.pdf',plot = pp5, width = 4, height = 4)
