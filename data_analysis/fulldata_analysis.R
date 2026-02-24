rm(list = ls())
library(Seurat) 
library(dplyr) 
library(ggplot2)
library(patchwork) 
library(ggpubr) 
library(RColorBrewer) 
library(grid)
library(ggrepel)
library(clusterProfiler)
# font_add("Arial", "C:/Windows/Fonts/arial.ttf")
# sysfonts::font_families()

# theme_set(
#   theme_minimal(base_family = "Arial", base_size = 18) +
#     theme(
#       text = element_text(face = "plain"),         
#       axis.title = element_text(face = "bold"),
#       axis.text = element_text(face = "plain"),
#       legend.title = element_text(face = "bold"),
#       legend.text = element_text(face = "plain"),
#       plot.title = element_text(face = "bold"),
#       strip.text = element_text(face = "plain")
#     )
# )
scPDAC <- readRDS("./Data/scPDAC_anno.rds")
#scPDAC_ref <- readRDS("./Data/scPDAC_anno_ref.rds")

scPDAC$celltype <- factor(scPDAC$celltype,levels = c("Acinar","Plasma","Neutrophils","Mast","Macrophage","Endothelial","Fibroblast","B cells","Ductal","NK cells","T cells"))
table(scPDAC$celltype)

# ====================
# ===== fig1.b =======
# ====================
library(tidydr)
library(ggsci)
library(paletteer)

scPDAC[['cellType']] = scPDAC@meta.data$c
umap = scPDAC@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(celltype = scPDAC@meta.data$celltype)
celltypepos <- umap %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = mean(UMAP_1),
    umap_2 = mean(UMAP_2))

my_col <- c("#E26CA8", "#8D78B7", "#E5007E", "#6593CE", "#31B6A9",
            "#70C9EB", "#E95D69", "#FBC89D", "#EDE346", "#F1832B", "#90C278")

p1.b <- DimPlot(scPDAC, label = FALSE, group.by = "celltype")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  geom_label_repel(aes(x = umap_1,y = umap_2,label = celltype,color = celltype), 
                   fontface = "plain",
                   data = celltypepos,
                   box.padding = 0.5) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Fig1/fig1b_1.pdf',plot = p1.b,width = 6, height = 5)

# ====================
# ===== fig1.e =======
# ====================

group_col <- c("Sensitive" = "#E95D69", "Resist" = "#5CACEE")
p1.b_2 <- DimPlot(scPDAC, cols = group_col, group.by = "group")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14, face = 2,hjust = 0.02),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 24)) +
  labs(color = "Group", title = NULL) +  # 设置图例标题
  theme(plot.title = element_text(size = 20)) +
  scale_color_manual(values = group_col, labels = c("Sensitive", "Resistant"))  # 明确指定颜色和标签
ggsave('./Result/Fig1/fig1e.pdf',plot = p1.b_2,width = 8, height = 5)

# ====================
# ===== fig1.g =======
# ====================

Idents(scPDAC) <- "celltype"
p1.g <- FeaturePlot(scPDAC,features = 'deg_num', label = TRUE, label.size = 6)+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) + 
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c("#6593CE",
                                    "#70C9EB",
                                    "#FFBB78FF",
                                    "#F3D5D8",
                                    "#EDA7A7",
                                    "#FF9896FF")) +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 18, face = "plain"), 
        legend.title = element_text(size = 24),  
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 16, face = 2, hjust = 0.02)) +
  labs(title = "",
       color = "DEG number")  

ggsave('./Result/Fig1/fig1g.pdf',plot = p1.g, width = 9, height = 6)

# ====================
# ===== fig1.c =======
# ====================
features <- list("T cells" = c("CD3D", "CD3E", "CD3G", "IL7R"),
                 "NK cells" = c("KLRD1", "NKG7", "CCL5", "GNLY"), 
                 "Ductal" = c("TSPAN8", "TFF2",  "AGR2", "KRT19", "MUC1", "EPCAM", "FXYD3", "CEACAM6", "KRT18"), 
                 "B cells" = c("MS4A1", "TNFRSF13C", "LY9", "BANK1"), 
                 "Fibroblast" = c("COL1A2", "DCN", "COL1A1", "LUM", "ACTA2", "COL3A1"),
                 "Endothelial" = c("PLVAP", "PECAM1"),
                 "Macrophage" = c("CD68", "APOE", "CXCL2", "C1QA", "CXCL3", "AIF1", "CD14", "FCGR3A"),
                 "Mast" = c("TPSB2", "MS4A2"),
                 "Neutrophils" = c("S100A9", "S100A8"),
                 "Plasma" = c("MZB1", "CD38", "IGHA1", "FKBP11"),
                 "Acinar" = c("PRSS1", "PRSS2"))
scPDAC$celltype <- factor(scPDAC$celltype,levels = rev(c("T cells", "NK cells", "Ductal", "B cells", "Fibroblast", "Endothelial", "Macrophage", "Mast", "Neutrophils", "Plasma", "Acinar")))

DefaultAssay(scPDAC) <- "RNA"
p1.c <- DotPlot(scPDAC, features = features, group.by = "celltype") + 
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),  
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12), 
        axis.text.y = element_text(size = 14)) +  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#330066','#336699',"#F3D5D8","#FF9896FF")) + 
  labs(x=NULL,y=NULL)


ggsave('./Result/Fig1/fig1c.pdf', plot = p1.c, width = 20, height = 5)

# ====================
# ===== fig1.d =======
# ====================
library(ggpubr)

portion <- prop.table(table(scPDAC@meta.data$celltype, scPDAC@meta.data$orig.ident), margin = 2)
portion <- as.data.frame(portion)
colnames(portion) <- c("celltype","orig.ident","ratio")
portion <- merge(portion,unique(scPDAC@meta.data[,c("orig.ident",'group'),drop=F]))
colnames(portion)[4] <- 'group'
portion$ratio <- round(portion$ratio*100,2)
portion$group <- factor(portion$group,levels = c('Sensitive','Resist'))
portion$celltype <- factor(portion$celltype, levels = c("T cells", "NK cells", "B cells", "Macrophage", "Plasma", "Neutrophils", "Ductal", "Acinar", "Endothelial", "Fibroblast", "Mast"))

p1.d <- ggplot(portion,aes(celltype,ratio,fill=group))+ 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Sensitive" = "#E95D69","Resist" = "#5CACEE"),
                    labels = c("Sensitive" = "Sensitive", "Resist" = "Resistant"))+
  labs(x = "", y = "", fill = "Group") +
  theme_bw()+  
  theme(text=element_text(family="sans"),
        axis.text.y = element_text(size=20, color = "black"),
        axis.text.x = element_text(size=22, angle = 45, hjust = 1, vjust = 1, color = "black"),        
        axis.title.y = element_text(size = 20, color = "black"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))  
ggsave('./Result/Fig1/fig1d.pdf', plot = p1.d, width = 11, height = 5)


# ====================
# ===== fig1.f =======
# ====================
#保存每种celltype的umap1和umap2维度
library(writexl)
umap_coords <- data.frame(scPDAC@reductions[["umap"]]@cell.embeddings,
                          "celltype" = scPDAC@meta.data$celltype,
                          "group" = scPDAC@meta.data$group)
for(celltype in unique(scPDAC$celltype)){
  umap_celltype_sen <- umap_coords[umap_coords$celltype == celltype & umap_coords$group == "Sensitive", c("UMAP_1", "UMAP_2")]
  umap_celltype_res <- umap_coords[umap_coords$celltype == celltype & umap_coords$group == "Resist", c("UMAP_1", "UMAP_2")]
  dataframes_list <- list(Sheet1 = umap_celltype_sen, Sheet2 = umap_celltype_res)
  write_xlsx(dataframes_list, path = paste("./Matlab/umap_coords/", celltype, ".xlsx", sep = ""))
}

# ..... 在matlab中计算kavaslue：
#run kstest.m/ Result in "./Data/ksvalue.xlsx"

library(readxl)
library(ggplot2)
library(dplyr)


ksvalue <- read_excel("./Data/ksvalue.xlsx")
ksvalue <- ksvalue %>% arrange(desc(ksvalue))
ksvalue$`p-value` <- -log10(ksvalue$`p-value`)
colnames(ksvalue)[which(colnames(ksvalue) == "p-value")] <- "-log10(p)"

p1.f_dot <- ggplot(ksvalue, aes(x = 0, y = reorder(celltype, ksvalue), size = ksvalue , fill = `-log10(p)`)) +
  geom_point(alpha = 0.7, shape = 21, stroke = 0.2) +
  scale_fill_gradientn( values = seq(0, 1, 0.25), colours = c(
    rgb(002, 038, 062, maxColorValue = 255),
    rgb(115, 186, 214, maxColorValue = 255),
    rgb(239, 065, 067, maxColorValue = 255),
    rgb(191, 030, 046, maxColorValue = 255))) +  
  labs(y = "celltype") +
  xlim(-0.1, 0.1) +  # 设置横轴显示范围
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12))
ggsave('./Result/Fig1/fig1f_dot.pdf',plot = p1.f_dot, width = 2.8, height = 6)

p1.f_dot <- ggplot(ksvalue, aes(x = reorder(celltype, ksvalue), y = 0, size = ksvalue, fill = `-log10(p)`)) +
  geom_point(alpha = 0.7, shape = 21, stroke = 0.2) +
  scale_fill_gradientn(values = seq(0, 1, 0.25), colours = c(
    rgb(2, 38, 62, maxColorValue = 255),
    rgb(115, 186, 214, maxColorValue = 255),
    rgb(239, 65, 67, maxColorValue = 255),
    rgb(191, 30, 46, maxColorValue = 255))) +
  labs(x = "celltype") +   # 横轴显示细胞类型
  ylim(-0.1, 0.1) +       # 设置纵轴显示范围
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.x = element_text(size = 12))

ggsave('./Result/Fig1/fig1f_dot_horizontal.pdf', plot = p1.f_dot, width = 6, height = 3)



