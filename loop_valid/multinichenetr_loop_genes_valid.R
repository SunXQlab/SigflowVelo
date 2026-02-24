rm(list = ls())
library(Seurat) 
library(dplyr) # %>%
library(ggplot2)
library(patchwork) 
library(ggpubr) # ggboxplot
library(RColorBrewer) # brewer.ylorrd
library(grid)
library(ggrepel)
library(clusterProfiler)
library(tidydr)
library(ggsci)
library(paletteer)

theme_set(theme_gray(base_family = "Arial"))

scPDAC_CCC <- readRDS("./Data/scPDAC_CCC.rds")
scPDAC_CCC$celltype <- as.character(scPDAC_CCC$celltype)
scPDAC_CCC$celltype[scPDAC_CCC$celltype == "Dutcal"] <- "Ductal"
scPDAC_CCC$celltype <- factor(scPDAC_CCC$celltype, levels = c("Macrophage", "T cells", "Ductal"))

## UMAP 
scPDAC_CCC[['cellType']] = scPDAC_CCC@meta.data$c
umap = scPDAC_CCC@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(celltype = scPDAC_CCC@meta.data$celltype)

celltypepos <- umap %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = mean(UMAP_1),
    umap_2 = mean(UMAP_2))

my_col <- c("#E95D69", "#6593CE", "#90C278")

p1.b <- DimPlot(scPDAC_CCC, label = FALSE, group.by = "celltype")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = celltype,color = celltype), 
                   fontface = "bold",
                   data = celltypepos,
                   box.padding = 0.5) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/fig7/umap.pdf',plot = p1.b,width = 6, height = 5)

## Degenes Dotplot
Idents(scPDAC_CCC) <- "celltype"

deglist <- list()
for(c in unique(scPDAC_CCC$celltype)){
  sce1 <- subset(scPDAC_CCC, celltype==c)
  degs <- FindMarkers(sce1,
                      ident.1 = 'Resist',
                      ident.2 = 'Sensitive',
                      assay = "RNA",
                      min.pct = 0.01,
                      logfc.threshold = 0,
                      group.by = 'group')
  degs <- cbind(gene=rownames(degs),degs)
  degs$celltype <- c
  deglist[[c]] <- degs
}
degs <- Reduce(rbind,deglist)
write.csv(degs,'./Data/celltype_group_degs.csv',quote = F,row.names = F)

degs <- read.csv("./Data/celltype_group_degs.csv")
degs$celltype[degs$celltype == "Dutcal"] <- "Ductal"

loop_degenes <- degs %>%
  filter(
    (celltype == "T cells" & gene %in% c("TNF", "SPP1", "PSAP", "SPP1", "GRN", "TNFRSF9", "TNFRSF1A", "TNFRSF1B", "TGFBR2", "ITGB1", "ITGA4", "CD74", "CD44", "IL1R2", "CCR5")) |
    (celltype == "Ductal" & gene %in% c("CD44", "SORT1", "ITGAV", "ITGB1", "ITGB5", "ITGB6", "ITGA5", "TGFBR1", "TGFBR2", "ITGAV", "ITGB1", "ITGB5", "ITGB6", "ITGA5", "CD44", "EGFR", "CD74", "ERBB2", "ERBB3")) |
    (celltype == "Macrophage" & gene %in% c("TGFB1", "SPP1", "MIF", "AREG", "TNFRSF1B", "CD44", "TNF", "TNFSF9", "TGFB1", "SPP1", "MIF", "IL1B", "CCL4"))
  )

loop_degenes <- loop_degenes %>%
  mutate(log10_pval = -log10(p_val))

loop_degenes$gene <- factor(loop_degenes$gene, levels = unique(loop_degenes$gene))
loop_degenes$celltype <- factor(loop_degenes$celltype, levels = unique(loop_degenes$celltype))

p3 <- ggplot(loop_degenes, aes(x = gene, y = celltype)) +
  geom_point(aes(size = log10_pval, color = avg_log2FC)) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(
    title = "",
    x = "",
    y = "",
    color = "avg_log2FC",
    size = "-log10(p_val)"
  ) +
  theme_bw(base_size = 13) + #coord_flip() +
  theme(panel.grid = element_blank(),  
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),  
        axis.text.y = element_text(size = 14)) + 
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#336699', "#FFC1C1","#FF9896FF",  "#E95D69")) 

ggsave('./Result/fig7/degenes_dotplot.pdf',plot = p3,width = 14, height = 3)

## FeaturePlot
p.2 <- FeaturePlot(scPDAC_CCC, features = c("TNF", "SPP1"), split.by = "group", max.cutoff = 3, 
            cols = c("#DCDCDC", "#E95D69", "#E95D69", "#E95D69"))

ggsave('./Result/fig7/features_plot.pdf',plot = p.2,width = 8, height = 6)
