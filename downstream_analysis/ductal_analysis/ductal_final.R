library(Seurat)
library(paletteer)
library(ggplot2)
library(dplyr)
library(tidydr)
library(ggrepel)
library(ggsci)
library(RColorBrewer)

scPDAC_ductal <- readRDS('./Data/ductal/scPDAC_ductal.rds')

DefaultAssay(scPDAC_dutcal)<-"RNA"
scPDAC_ductal_normalized <- NormalizeData(scPDAC_ductal, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

saveRDS(scPDAC_ductal_normalized, file = "./Data/ductal/scPDAC_ductal_normalized.rds")

#========== umap group ==============#
my_col <- c("#E95D69", "#5CACEE" )

pp4 <- DimPlot(scPDAC_ductal, label = F, group.by = "group")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Ductal/umap_group.pdf',plot = pp4, width = 6, height = 6)

#============ umap anno ==============#
my_col <- c("#E26CA8", "#8D78B7", "#E5007E", "#6593CE", "#31B6A9",
            "#70C9EB", "#E95D69", "#FBC89D", "#EDE346")
pp4 <- DimPlot(scPDAC_ductal, label = T, group.by = "anno", label.size = 5)+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=my_col)
ggsave('./Result/Ductal/umap_anno.pdf',plot = pp4, width = 6, height = 6)


# =========== E/M score ==============#
genesets <- list()
genesets$E <- c("OCLN", "DSP", "ESRP1", "ESRP2", 
                "MUC1", "TJP3", "LLGL2", 
                "ST14", "GRHL2", "OVOL2", "DSC2")

genesets$M <- c("VIM", "FN1", "TWIST1", "ZEB2", "FOXC2", 
                "MMP2", "MMP9", "FGF2", "PDGFRB", "ACTA2", "COL1A1", "COL3A1", "COL5A1", 
                "SPARC", "ITGB1", "TNC", "POSTN", "SERPINE1")



scPDAC_ductal <- AddModuleScore(scPDAC_ductal,
                                assay = 'RNA',
                                features = genesets,
                                ctrl = 5)
name_index <- ncol(scPDAC_ductal@meta.data)
colnames(scPDAC_ductal@meta.data)[(name_index-1):name_index] <- c('E','M')

# EMT Score
E_genes <-  c("OCLN", "DSP", "ESRP1", "ESRP2", 
              "MUC1", "TJP3", "LLGL2", 
              "ST14", "GRHL2", "OVOL2", "DSC2")

M_genes <-  c("VIM", "FN1", "TWIST1", "ZEB2", "FOXC2", 
              "MMP2", "MMP9", "FGF2", "PDGFRB", "ACTA2", "COL1A1", "COL3A1", "COL5A1", 
              "SPARC", "ITGB1", "TNC", "POSTN", "FSP1", "SERPINE1")
scPDAC_ductal <- AddModuleScore(scPDAC_ductal, features = list(E_genes), name = "E_Score")
scPDAC_ductal <- AddModuleScore(scPDAC_ductal, features = list(M_genes), name = "M_Score")

scPDAC_ductal$EMT_Score <- scPDAC_ductal$M_Score1 - scPDAC_ductal$E_Score1

# 计算p-value
E_pvalue <- wilcox.test(E ~ group, data = scPDAC_ductal@meta.data)$p.value
print(E_pvalue)
M_pvalue <- wilcox.test(M ~ group, data = scPDAC_ductal@meta.data)$p.value
print(M_pvalue)
EMT_pvalue <- wilcox.test(EMT_Score ~ group, data = scPDAC_ductal@meta.data)$p.value
print(EMT_pvalue)

group_levels <- levels(scPDAC_ductal@meta.data$group)
x_mid <- 1.5 # 分组的中间位置
y_max <- 1.2 # 设置 y 轴上方的空白区域

p1 <- ggplot(scPDAC_ductal@meta.data, aes(x = group, y = E, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  # 使用小提琴图
  geom_boxplot(width = 0.1,          # 箱线图宽度
               fill = "white",       
               alpha = 0.8,         
               outlier.shape = NA) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-1, 1.4)) +  # 只限制y轴
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) +  
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  xlab('E score') +
  theme(legend.position = "none") +
  annotate("text", x = x_mid, y = 1.4, label = paste0("p = ", signif(E_pvalue, digits = 3)), size = 5, hjust = 0.5)

p2 <- ggplot(scPDAC_ductal@meta.data, aes(x = group, y = M, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  # 使用小提琴图
  geom_boxplot(width = 0.1,          # 箱线图宽度
               fill = "white",       
               alpha = 0.8,         
               outlier.shape = NA) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-0.5, 1.2)) +  # 只限制y轴（Cytotoxic分数）
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) + 
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  xlab('M score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  theme(legend.position = "none") +
  annotate("text", x = x_mid, y = y_max, label = paste0("p = ", signif(M_pvalue, digits = 3)), size = 5, hjust = 0.5)

p3 <- ggplot(scPDAC_ductal@meta.data, aes(x = group, y = EMT_Score, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  # 使用小提琴图
  geom_boxplot(aes(fill = group), width = 0.1, fill = "white", alpha = 0.8, outlier.shape = NA) +  # 箱线图不映射到图例
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  # 只限制y轴（Cytotoxic分数）
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) + 
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  xlab('EMT score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  annotate("text", x = x_mid, y = 2, label = paste0("p = ", signif(EMT_pvalue, digits = 3)), size = 5, hjust = 0.5)


p <- p1 + p2 + p3
ggsave('./Result/Ductal/EM_score_violin.pdf', plot = p, width = 12, height = 2.67, dpi = 300)

# EMT Score随Pseudotime变化
E_genes <-  c("CLDN1", "OCLN", "DSP", "ESRP1", "ESRP2", 
              "MUC1", "TJP3", "LLGL2", 
              "ST14", "GRHL2", "OVOL2", "DSC2", "DSG2")

M_genes <-  c("VIM", "FN1", "TWIST1", "ZEB2", "FOXC2", 
              "MMP2", "MMP9", "FGF2", "PDGFRB", "ACTA2", "COL1A1", "COL3A1", "COL5A1", 
              "SPARC", "ITGB1", "TNC", "POSTN", "FSP1", "SERPINE1")

#============ 置信区间：Abnormal和Ressit的关系 ==============#
library(Hmisc)
library(dplyr)
library(ggplot2)

cluster_group <- scPDAC_ductal@meta.data[, c('seurat_clusters', 'group')]

count_group <- cluster_group %>%
  group_by(seurat_clusters, group) %>%
  summarise(count = n()) %>%
  ungroup()

proportion_group <- count_group %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

new_row <- data.frame(
  seurat_clusters = factor(6),
  group = "Resist",
  count = 0,
  proportion = 0
)

proportion_group <- rbind(proportion_group, new_row)

resist_data <- proportion_group %>%
  filter(group == "Resist") %>%
  arrange(desc(proportion)) %>%
  rename(resist_proportion = proportion)

cluster_isNormal <- scPDAC_ductal@meta.data[, c('seurat_clusters', 'isNormal')]

count_isNormal <- cluster_isNormal %>%
  group_by(seurat_clusters, isNormal) %>%
  summarise(count = n()) %>%
  ungroup()

proportion_isNormal <- count_isNormal %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

abnormal_data <- proportion_isNormal %>%
  filter(isNormal == "Abnormal") %>%
  arrange(desc(proportion)) %>%
  rename(abnormal_proportion = proportion)

merged_data <- left_join(resist_data, abnormal_data, by = "seurat_clusters")
merged_data <- merged_data %>%
  select(-count.x, -count.y, -group, -isNormal) %>%
  arrange(desc(resist_proportion))

cor_test_result <- cor.test(
  x = merged_data$resist_proportion,
  y = merged_data$abnormal_proportion,
  method = "pearson"
)

pearson_cor <- cor_test_result$estimate
p_value <- cor_test_result$p.value

cor_label <- paste("Pearson r =", round(pearson_cor, 2))

print(paste("Pearson r =", round(pearson_cor, 4)))
print(paste("p-value =", signif(p_value, 3)))

# 绘制散点图，并添加带状区域（置信区间）和 相关系数
p2.e <- ggplot(merged_data, aes(x = resist_proportion, y = abnormal_proportion, color = seurat_clusters)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "#E95D69", fill = alpha("#AEEEEE", 0.2)) +  # 添加带状区域
  geom_text(aes(label = seurat_clusters), vjust = -1) +
  annotate("text", x = max(merged_data$resist_proportion) * 0.8, y = max(merged_data$abnormal_proportion) * 0.8, 
           label = cor_label, color = "#E95D69", size = 5) +  # 添加 R-square
  labs(title = "",
       x = "Proportion of Resistant Cells",
       y = "Proportion of Malignant Cells ") +
  guides(color = FALSE) +  # 不显示颜色图例
  theme_minimal() +
  scale_color_manual(values = my_col) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"))  # 显示坐标轴，设置坐标轴样式
ggsave('./Result/Ductal/resist_vs_abnormal.pdf', plot = p2.e, width = 5, height = 3.5)


#============ CNV Vliont Plot ==============#
library(rjags)
library(infercnv)
infercnv_obj <- readRDS('./Data/ductal/run.final.infercnv_obj')

cnv_mat <- as.data.frame(t(infercnv_obj@expr.data))

cnvMean <- function(data) {
  data <- data %>% as.matrix() 
  abs_diff <- abs(1 - data)
  cnv_mean <- rowMeans(abs_diff)
  return(as.data.frame(cnv_mean))
}

cnv_mean <- cnvMean(cnv_mat)

scPDAC <- readRDS("./Data/scPDAC_anno.rds")
celltype_info <- scPDAC@meta.data$celltype
names(celltype_info) <- rownames(scPDAC@meta.data)

common_cells <- intersect(rownames(cnv_mean), names(celltype_info))

cnv_mean <- cnv_mean[common_cells, , drop = FALSE]  # 保持为 dataframe

cnv_mean$celltype <- celltype_info[common_cells]

cnv_mean$group <- ifelse(cnv_mean$celltype == "Ductal", 
                         "Ductal", 
                         "Other Celltypes")

# Calculate p-value using Wilcoxon test
p_value <- wilcox.test(cnv_mean ~ group, data = cnv_mean)$p.value

# Create the plot
p <- ggplot(cnv_mean, aes(x = group, y = cnv_mean, fill = group)) +
  # Violin plot with trimmed edges and transparency
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  
  # Boxplot inside violin with custom settings
  geom_boxplot(width = 0.15, 
               fill = "white", 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  
  # Custom color scheme
  scale_fill_manual(values = c("Ductal" = "#E5007E", "Other Celltypes" = "#70C9EB")) +
  scale_colour_manual(values = c("Ductal" = "#E5007E", "Other Celltypes" = "#70C9EB")) + 
  
  # Axis labels and title
  labs(
    title = "",
    x = "",
    y = "CNV Score"
  ) +
  
  # Theme customization
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    panel.grid = element_blank(),  # Remove all grid lines
    panel.background = element_blank(),  # Clear background
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  ) +
  
  # p-value annotation
  annotate("text", 
           x = 1.5, 
           y = max(cnv_mean$cnv_mean) * 0.95,
           label = paste0("p = ", format.pval(p_value, digits = 2)),
           size = 4.5)

# Save the plot
ggsave("./Result/Ductal/cnv_violin_Ductal_vs_Others.pdf", 
       plot = p, 
       width = 5, 
       height = 6)


## 查看最右下角的cell id

umap_coords <- Embeddings(scPDAC_ductal, "umap")

# 转为数据框
umap_df <- as.data.frame(umap_coords)

# 计算“右下角”分数：右（UMAP_1大）- 下（UMAP_2小）
umap_df$score <- umap_df$UMAP_1 - umap_df$UMAP_2

# 找到得分最大的细胞（即最右下）
bottom_right_cell <- rownames(umap_df)[which.max(umap_df$score)]
bottom_right_cell

library(ggplot2)

umap_df <- as.data.frame(umap_coords)
umap_df$cell_id <- rownames(umap_df)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_point(data = subset(umap_df, cell_id == bottom_right_cell),
             color = "red", size = 2) +
  theme_bw()
