rm(list = ls())
library(Seurat)
library(paletteer)
library(ggplot2)
library(dplyr)

## load dataT
scPDAC_T <- readRDS("./Data/Tcells/scPDAC_T.rds")

# Dotplot
features <- c(#"CCL4L2", "CCL4", "IFNG",  # CD8(IFNG)
  "ANXA1", "HSPD1", "HSP90AA1", "MYADM", # CD4(CCL20)
  "FOS", "NFKBIA", # CD4(FOS)
  #"TUBA4A",  # CD8(GZMK)
  #"RPS8", "RPS3A", "RPS6","CCR7", # CD4(CCR7) TN
  "RPS3A", "CCR7", # CD4(CCR7) TN
  #"FTH1",  # CD8(FTH1)
  #"TRAV1-2", "TNF", "NCR3",  # CD8 KLRB1+
  "BATF", "TNFRSF4", "PMAIP1", "IL2RA", "TIGIT", "RGS1", "FOXP3", # CD4 T_reg
  "MT1X", "MT2A","GZMA", # CD8(GZMA)
  #"CD52", "LTB", "IL7R", # CD4(LTB)
  #"KIAA1551", "GZMK", "EOMES"# CD8(EOMES)
  "STMN1", "KIAA1551" # CD8 proliferating
)

DefaultAssay(scPDAC_T)<-"RNA"
scPDAC_T$celltype <- factor(scPDAC_T$celltype, levels = c("CD4 CCL20+", "CD4 FOS", "CD4 Naive", "CD4 Treg", "CD8 GZMA+", "CD8 proliferating"))

p <- DotPlot(scPDAC_T, features = features, group.by = "celltype") + 
  theme_bw(base_size = 13) + 
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5, size = 14),
        axis.text.y = element_text(size = 14)) + 
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c("white", alpha("#BBFFFF",0.2), alpha("#AEEEEE", 0.2), alpha("#98F5FF", 0.2),"#F3D5D8", "#FF9896FF","#FF6A6A")) + 
  labs(x=NULL,y=NULL) 
ggsave("./Result/Tcells/Dotplot.pdf", plot = p,width = 12,height = 4)

scPDAC_T$seurat_clusters <- factor(scPDAC_T$seurat_clusters, levels = c("0", "4", "8", "3", "7", "1", "2", "5", "9", "6"))

#markers
scPDAC_T.markers <- FindAllMarkers(scPDAC_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(scPDAC_T.markers, file = "./Data/Tcells/T_allmarkers.rds")

#============ Fig2a ==============#
scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')
scPDAC_T.markers <- readRDS('./Data/Tcells/T_allmarkers.rds')
DefaultAssay(scPDAC_T) <- "RNA"

library(tidydr)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(paletteer)

#========== group ==============#
my_col <- c("#E95D69", "#5CACEE" )

pp4 <- DimPlot(scPDAC_T, label = F, group.by = "group")+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend()+
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Tcells/umap_group_nolegend.pdf',plot = pp4, width = 6, height = 6)

#============ celltype ==============#
my_col <- c("#6593CE","#70C9EB", "#EF8B67", "#F5EBAE", "#F7FBC9", "#F0C284")
pp4 <- DimPlot(scPDAC_T, label = T, group.by = "celltype", label.size = 6)+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Tcells/umap_celltype.pdf',plot = pp4, width = 6, height = 6)


#============ CD4&CD8 ==============#
my_col <- c("#7B68EE", "#FF8247")
pp5 <- DimPlot(scPDAC_T, label = T, group.by = "subtype", label.size = 7)+ 
  theme_dr(xlength = 0.2, 
           ylength = 0.25,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) +
  NoLegend() +
  labs(title = NULL) +
  scale_color_manual(values=rev(my_col))
ggsave('./Result/Tcells/umap_subtype.pdf',plot = pp5, width = 6, height = 6)

# ===== group cell counts ratio in different clusters ===== #
cluster_group <- scPDAC_T@meta.data[,c('seurat_clusters','group')]

count_data <- cluster_group %>%
  group_by(seurat_clusters, group) %>%
  summarise(count = n()) %>%
  ungroup()

proportion_data <- count_data %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

colors <- c("Resist" = "#5CACEE", "Sensitive" = "#E95D69")

pp3.pie <- ggplot(proportion_data, aes(x = "", y = proportion, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + 
  labs(fill = "group") +
  ggtitle("Proportion of Groups in Seurat Clusters") +
  facet_wrap(~seurat_clusters) +
  scale_fill_manual(values = colors)

ggsave('./Result/Tcells/cluster_group_pie.pdf',plot = pp3.pie, height = 6, width = 6)

#=================================#
#============ Fig2b ==============#
#=================================#
library(ggplotify)
source("pheatmap_add_gene_name.R")

all_genes <- rownames(scPDAC@assays$RNA)
#Growth Factor
growthfactor <- c(grep("^TGF", all_genes, value = TRUE), 
                  grep("^BMP", all_genes, value = TRUE),
                  grep("^PDGF", all_genes, value = TRUE),
                  grep("^FGF", all_genes, value = TRUE),
                  grep("^EGF", all_genes, value = TRUE), 
                  grep("^VEGf", all_genes, value = TRUE),
                  grep("^NGF", all_genes, value = TRUE),
                  grep("^IGF", all_genes, value = TRUE))

# ======= Group
DefaultAssay(scPDAC_T)<-"RNA"
scPDAC_T$newgroup <- paste(scPDAC_T$seurat_clusters, scPDAC_T$group,sep = '_')
bs <- split(colnames(scPDAC_T), scPDAC_T[['newgroup']])
exprSet <- do.call(
  cbind,lapply(names(bs), function(x){
    kp <- colnames(scPDAC_T) %in% bs[[x]]
    rowSums(as.matrix(scPDAC_T@assays$RNA@counts[, kp]))
  })
)
colnames(exprSet) <- names(bs)
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[which(rowSums(exprSet) > 0),]
# exprSet <- exprSet[rownames(exprSet) %in% sce.markers$gene,]
anno <- unique(scPDAC_T@meta.data[,c('newgroup','group','seurat_clusters')])
rownames(anno) <- anno$newgroup
anno$newgroup <- NULL
ann_colors = list(
  seurat_clusters = c("0" = "#1B9E77", "1" = "#D95F02", "2" ="#7570B3", "3" ="#E7298A", "4" ="#66A61E", "5" = "#E41A1C", "6" = "#377EB8", "7" = "#4DAF4A", "8" = "#F3C300"),
  group = c(Resist = "#5CACEE", Sensitive = "#E95D69")
)

p1 <- pheatmap(exprSet,
               cluster_cols = F,
               annotation_col = anno,
               treeheight_row = 0,
               annotation_colors = ann_colors,
               color = colorRampPalette(c('#084594','white','#C10525'))(100),
               show_colnames = F,
               scale = "row")
p2 <- add.flag(p1,
               kept.labels = growthfactor,
               repel.degree = 0.2)
p2 <- as.ggplot(p2)
ggsave('./Result/Tcells/geneExpression_growthfactor.pdf',plot = p2,width = 6.48,height = 10)


# ======= Patient Sensitive
scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')
DefaultAssay(scPDAC_T)<-"RNA"
target_orig_ident <- c("P1-T", "P2-T", "P3-T", "P5-T", "P6-T")
scPDAC_T <- subset(scPDAC_T, orig.ident %in% target_orig_ident)
scPDAC_T$newgroup <- paste(scPDAC_T$seurat_clusters, scPDAC_T$orig.ident, sep = '_')

bs <- split(colnames(scPDAC_T), scPDAC_T[['newgroup']])

exprSet <- do.call(
  cbind, lapply(names(bs), function(x) {
    kp <- colnames(scPDAC_T) %in% bs[[x]]
    rowSums(as.matrix(scPDAC_T@assays$RNA@counts[, kp]))
  })
)

colnames(exprSet) <- names(bs)

exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[which(rowSums(exprSet) > 0),]

anno <- unique(scPDAC_T@meta.data[, c('newgroup', 'orig.ident', 'seurat_clusters')])
rownames(anno) <- anno$newgroup
anno$newgroup <- NULL

ann_colors = list(
  seurat_clusters = c("0" = "#1B9E77", "1" = "#D95F02", "2" = "#7570B3", "3" = "#E7298A", "4" = "#66A61E", "5" = "#E41A1C", "6" = "#377EB8", "7" = "#4DAF4A", "8" = "#F3C300"),
  orig.ident = c("P1-T" = "#5CACEE", "P2-T" = "#E95D69", "P3-T" = "#FFD700", "P5-T" = "#ADFF2F", "P6-T" = "#FF69B4")
)

p1 <- pheatmap(exprSet,
               cluster_cols = F,
               annotation_col = anno,
               treeheight_row = 0,
               annotation_colors = ann_colors,
               color = colorRampPalette(c('#084594','white','#C10525'))(100),
               show_colnames = F,
               scale = "row")
p2 <- add.flag(p1,
               kept.labels = growthfactor,
               repel.degree = 0.2)
p2 <- as.ggplot(p2)
ggsave('./Result/Tcells/geneExpression_Sen.pdf',plot = p2,width = 6.48,height = 8)

# ======= Patient Resist
scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')
DefaultAssay(scPDAC_T)<-"RNA"
target_orig_ident <- c("P7-T", "P8-T")
scPDAC_T <- subset(scPDAC_T, orig.ident %in% target_orig_ident)
scPDAC_T$newgroup <- paste(scPDAC_T$seurat_clusters, scPDAC_T$orig.ident, sep = '_')

bs <- split(colnames(scPDAC_T), scPDAC_T[['newgroup']])

exprSet <- do.call(
  cbind, lapply(names(bs), function(x) {
    kp <- colnames(scPDAC_T) %in% bs[[x]]
    rowSums(as.matrix(scPDAC_T@assays$RNA@counts[, kp]))
  })
)

colnames(exprSet) <- names(bs)

exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[which(rowSums(exprSet) > 0),]

anno <- unique(scPDAC_T@meta.data[, c('newgroup', 'orig.ident', 'seurat_clusters')])
rownames(anno) <- anno$newgroup
anno$newgroup <- NULL

ann_colors = list(
  seurat_clusters = c("0" = "#1B9E77", "1" = "#D95F02", "2" = "#7570B3", "3" = "#E7298A", "4" = "#66A61E", "5" = "#E41A1C", "6" = "#377EB8", "7" = "#4DAF4A", "8" = "#F3C300"),
  orig.ident = c("P7-T" = "#FFD700", "P8-T" = "#FF69B4")
)

p1 <- pheatmap(exprSet,
               cluster_cols = F,
               annotation_col = anno,
               treeheight_row = 0,
               annotation_colors = ann_colors,
               color = colorRampPalette(c('#084594','white','#C10525'))(100),
               show_colnames = F,
               scale = "row")
p2 <- add.flag(p1,
               kept.labels = growthfactor,
               repel.degree = 0.2)
p2 <- as.ggplot(p2)
ggsave('./Result/Tcells/geneExpression_Res.pdf',plot = p2,width = 6.48,height = 8)

#=================================#
#========== Annoate ==============#
#=================================#
library(dplyr)
# scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')
# scPDAC_T.markers <- readRDS('./Data/Tcells/T_allmarkers.rds')


top_markers <- scPDAC_T.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC)

DefaultAssay(scPDAC_T)<-"RNA"
anno <- data.frame(  
  seurat_clusters = scPDAC_T@meta.data$seurat_clusters,  
  celltype = NA 
)  
anno$celltype[anno$seurat_clusters %in% c(0)] <- "CCL5+ Tcell"  
anno$celltype[anno$seurat_clusters %in% c(1)] <- "	TNFRSF4+ Tcell"  
anno$celltype[anno$seurat_clusters %in% c(2)] <- "HIST1H1D+ Tcell"  
anno$celltype[anno$seurat_clusters %in% c(3)] <- "	IGLV2-14+ Tcell"  
anno$celltype[anno$seurat_clusters %in% c(4)] <- "CD8B+ Tcell"  

scPDAC_T@meta.data$anno <- anno$celltype
p2a_anno <- DimPlot(scPDAC_T,group.by = 'anno', label = TRUE) + NoLegend()
ggsave('./Result/Tcells/umap_anno.pdf',plot = p2a_anno,width = 6,height = 6)


#====================================================#
#============ growthfactor expression ==============#
#====================================================#
library(ggpubr)
scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')
#growthfactor
all_genes <- rownames(scPDAC@assays$RNA)
#趋化因子
growthfactor <- c(grep("^TGF", all_genes, value = TRUE), 
                  grep("^BMP", all_genes, value = TRUE),
                  grep("^PDGF", all_genes, value = TRUE),
                  grep("^FGF", all_genes, value = TRUE),
                  grep("^EGF", all_genes, value = TRUE), 
                  grep("^VEGf", all_genes, value = TRUE),
                  grep("^NGF", all_genes, value = TRUE),
                  grep("^IGF", all_genes, value = TRUE))

Idents(scPDAC_T) <- scPDAC_T$group
#计算基因的平均表达量
groupExp <- AverageExpression(scPDAC_T)
groupExp <- groupExp[["RNA"]]
#找到groupExp中存在的基因
gMeanExp <- groupExp[growthfactor, ]
gMeanExp <- t(apply(gMeanExp, 1, function(x){(x - min(x)) / max(x)}))
gMeanExp <- reshape2::melt(t(gMeanExp)) %>% 
  dplyr::rename(Sample="Var1", Gene="Var2")
gMeanExp$Sample <- factor(gMeanExp$Sample)

cols <- setNames(c("#E95D69", "#5CACEE"),
                 c("Sensitive", "Resist"))

p <- ggboxplot(gMeanExp, x = "Sample", y = "value",color = "Sample", orientation = "horizontal", add = "jitter", ylab = "Scaled mean", palette = cols) + ggtitle("growthFactor")
ggsave('./Result/Tcells/growthFactorExpression.pdf',plot = p,width = 5,height = 5)

#=================================#
#============ GSVA ==============#
#=================================#
library(msigdbr)
library(tidyverse)
library(GSVA)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

scPDAC_T <- readRDS('./Data/Tcells/scPDAC_T.rds')

#============ 
DefaultAssay(scPDAC_T) <- "RNA"
Idents(scPDAC_T) <- 'group'
bs <- split(colnames(scPDAC_T),scPDAC_T[['group']])
exprSet <- do.call(
  cbind,lapply(names(bs), function(x){
    kp <- colnames(scPDAC_T) %in% bs[[x]]
    rowSums(as.matrix(scPDAC_T@assays$RNA@data[, kp]))
  })
)
colnames(exprSet) <- names(bs)
exprSet <- as.data.frame(exprSet)
df <- exprSet[which(rowSums(exprSet) > 0),]
h <- msigdbr(species = "Homo sapiens", category = "H")
h$gs_name <- gsub('HALLMARK_','',h$gs_name)
# h <- h %>% filter(str_detect(gs_name,'_DN|_UP'))
h <- dplyr::select(h, gs_name, gene_symbol) %>% 
  as.data.frame %>%
  split(., .$gs_name) %>%
  lapply(., function(x)(x$gene_symbol))
gsvas <- gsva(gsvaParam(as.matrix(df), h))
gsvas <- as.data.frame(gsvas)
write.csv(gsvas,'./Result/Tcells/hallmark_gsva_score.csv',quote = F)
pdf('./Result/Tcells/gsva.pdf',height = 10,width = 8)
draw(
  Heatmap(
    gsvas,
    name = "GSVA Scores",
    col  = colorRamp2(seq(from=-0.5,to=0.3,length=11),rev(brewer.pal(11, "Spectral"))),
    use_raster = T,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 12),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      legend_width = unit(6, "cm")),
    cluster_rows = TRUE, #
    show_row_dend = F,
    cluster_columns = FALSE),
  heatmap_legend_side = "top"
)
dev.off()


library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
load('D:/Bio/PDAC/PDAC/Result/Tcells/pt.matrix.rda')

pdf('D:/Bio/PDAC/PDAC/Result/Tcells/tmp.pdf')
Heatmap(
  pt.matrix,
  name = "z-score",
  col  = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  use_raster = T,
  show_row_names  = FALSE, #
  show_column_names  = FALSE,
  row_names_gp = gpar(fontsize = 10),
  # clustering_distance_rows = 'manhattan',
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  row_split = 3,
  row_title = c('M1', 'M2', 'M3'),
  cluster_rows = TRUE,
  show_row_dend = F,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE)
dev.off()

p2 <- Heatmap(
  pt.matrix,
  name = "z-score",
  col  = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  use_raster = T,
  show_row_names  = FALSE, #
  show_column_names  = FALSE,
  row_names_gp = gpar(fontsize = 10),
  # clustering_distance_rows = 'manhattan',
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  row_split = 3,
  row_title = c('M1', 'M2', 'M3'),
  cluster_rows = TRUE,
  show_row_dend = F,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE)
clusterlist = row_order(p2)
clu_df <- lapply(1:length(clusterlist), function(i){
  out <- data.frame(GeneID = rownames(pt.matrix)[clusterlist[[i]]],
                    Cluster = paste0("M", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  do.call(rbind, .)
write.csv(clu_df,'D:/Bio/PDAC/PDAC/Result/Tcells/monocle_pheatmap_cut_genes.csv',quote = F,row.names = F)

clu_df <- read.csv('D:/Bio/PDAC/PDAC/Result/Tcells/monocle_pheatmap_cut_genes.csv', header = TRUE)
GOS <- NULL
for(m in unique(clu_df$Cluster)){
  eg <- bitr(clu_df[clu_df$Cluster==m,'GeneID'], fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(eg$ENTREZID, 
                 OrgDb = org.Hs.eg.db, 
                 ont='BP',
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 1, 
                 qvalueCutoff = 1,
                 keyType = 'ENTREZID')
  go <- setReadable(go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  deggo <- go@result
  bgcount <- as.numeric(unlist(strsplit(go@result[1,3],'/'))[2])
  deggo$GeneRatios <- deggo$Count/bgcount
  deggo$cluster <- m
  GOS <- rbind(GOS,deggo)
}
GOS <- GOS[GOS$p.adjust<=0.05 & GOS$GeneRatios>=0.05,]
write.csv(GOS,'D:/Bio/PDAC/PDAC/Result/Tcells/monocle_pheatmap_cut_genes_GO.csv',quote = F)

GOS <- read.csv('D:/Bio/PDAC/PDAC/Result/Tcells/monocle_pheatmap_cut_genes_GO.csv', header = TRUE)
pdf('D:/Bio/PDAC/PDAC/Result/Tcells/monocle_pheatmap_GO.pdf',width = 5.56,height = 7.1)
Heatmap(
  pt.matrix,
  name = "z-score",
  col  = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  use_raster = T,
  show_row_names  = FALSE, #
  show_column_names  = FALSE,
  row_names_gp = gpar(fontsize = 10),
  # clustering_distance_rows = 'manhattan',
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  row_split = 3,
  row_title = c('regulation of T cell activation',
                'regulation of translation', 'RNA splicing'), #根据富集结果选择展示
  cluster_rows = TRUE,
  show_row_dend = F,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE)
dev.off()

go_gene_mapping <- data.frame(GO_Term = character(), Gene = character(), Description = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(GOS)) {
  go_term <- GOS$ID[i]  
  description <- GOS$Description[i] 
  genes_in_go <- GOS$geneID[i]  
  gene_list <- unlist(strsplit(genes_in_go, "/"))  
  
  for (gene in gene_list) {
    go_gene_mapping <- rbind(go_gene_mapping, data.frame(GO_Term = go_term, Gene = gene, Description = description, stringsAsFactors = FALSE))
  }
}

write.csv(go_gene_mapping, 'D:/Bio/PDAC/PDAC/Result/Tcells/GO_Gene_Enrichment.csv', row.names = FALSE, quote = FALSE)


## Exhausted Score
library(Seurat)
library(paletteer)
library(ggplot2)
library(patchwork)

scPDAC_T <- readRDS("./Data/Tcells/scPDAC_T.rds")
# T exhausted & cytotoxic score
genesets <- list()
genesets$exhausted <- c('PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2', 'BTLA', 'ENO1', 'GAPDH', 'PPAR', 'AMPK', 
                        'TOX', 'EOMES', 'TBX21', 'NR4A1', 'NR4A2', 'NR4A3', 'BATF',
                        'IFNG', 'IL2', 'CXCL13', 'CCR7', 'CD95', 'BCL2', 'BCL2L11', 'CDKN1A')
genesets$cytotoxic <- c('IL2', 'GZMA', 'GNLY', 'PRF1', 'GZMB', 'GZMK', 'IFNG', 'NKG7')
scPDAC_T <- AddModuleScore(scPDAC_T,
                           assay = 'RNA',
                           features = genesets,
                           ctrl = 5)
name_index <- ncol(scPDAC_T@meta.data)
colnames(scPDAC_T@meta.data)[(name_index-1):name_index] <- c('exhausted','cytotoxic')

# 计算p-value
exhausted_pvalue <- wilcox.test(exhausted ~ group, data = scPDAC_T@meta.data)$p.value
print(exhausted_pvalue)
cytotoxic_pvalue <- wilcox.test(cytotoxic ~ group, data = scPDAC_T@meta.data)$p.value
print(cytotoxic_pvalue)

group_levels <- levels(scPDAC_T@meta.data$group)
x_mid <- 1.5 
y_max <- 2 

p1 <- ggplot(scPDAC_T@meta.data, aes(x = group, y = exhausted, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-1, 2)) +  
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) +  
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  xlab('Exhausted score') +
  theme(legend.position = "none") +
  annotate("text", x = x_mid, y = y_max, label = paste0("p = ", signif(exhausted_pvalue, digits = 3)), size = 5, hjust = 0.5)+
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.1,
    color = "#2E8B57",
    size = 0.2
  )

p2 <- ggplot(scPDAC_T@meta.data, aes(x = group, y = cytotoxic, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-1, 2)) +  
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) + 
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  xlab('Cytotoxic score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  annotate("text", x = x_mid, y = y_max, label = paste0("p = ", signif(cytotoxic_pvalue, digits = 3)), size = 5, hjust = 0.5)+
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.1,
    color = "#2E8B57",
    size = 0.2
  )

p <- p1 + p2
ggsave('./Result/Tcells/Exhausted_Cytotoxic_score_violin_T.pdf', plot = p, width = 7.86, height = 2.67, dpi = 300)


p1 <- ggplot(scPDAC_T@meta.data, aes(x = celltype, y = exhausted, colour = celltype, fill = celltype)) +
  geom_violin(size = 1, alpha = 0.3) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) + 
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  xlab('Exhausted score') +
  theme(legend.position = "none")
p2 <- ggplot(scPDAC_T@meta.data, aes(x = celltype, y = cytotoxic, colour = celltype, fill = celltype)) +
  geom_violin(size = 1, alpha = 0.3) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  
  xlab('Cytotoxic score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14))
p <- p1 + p2
ggsave('./Result/Tcells/Exhausted_Cytotoxic_score_violin_T_celltype.tiff', plot = p, width = 15.86, height = 2.67, dpi = 300)


# CD4 naiveness & Treg score
CD4_T <- rownames(scPDAC_T@meta.data[scPDAC_T@meta.data$subtype %in% c("CD4"),])
scPDAC_T_CD4 <- subset(scPDAC_T,cells=CD4_T)
genesets <- list()
genesets$naiveness <- c('CCR7', 'CD62L', 'LEF1', 'TCF7','IL7R','CD27','FOXP1', 'IKZF1', 'ZNF683', 'BACH2', 'CD28', 'PTPN22', 'KLF2')
genesets$treg <- c('FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'LRRC32', 'TGFB1', 'IL10', 'IL35', 'ENTPD1', 'NT5E', 'TNFRSF18', 'CCR4', 'CCR8', 'STAT5', 'BACH2', 'PRDM1')
scPDAC_T_CD4 <- AddModuleScore(scPDAC_T_CD4,
                               assay = 'RNA',
                               features = genesets,
                               ctrl = 5)
name_index <- ncol(scPDAC_T_CD4@meta.data)
colnames(scPDAC_T_CD4@meta.data)[(name_index-1):name_index] <- c('naiveness','treg')

naiveness_pvalue <- wilcox.test(naiveness ~ group, data = scPDAC_T_CD4@meta.data)$p.value
print(naiveness_pvalue)
treg_pvalue <- wilcox.test(treg ~ group, data = scPDAC_T_CD4@meta.data)$p.value
print(treg_pvalue)

p1 <- ggplot(scPDAC_T_CD4@meta.data, aes(x = group, y = naiveness, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  # 使用小提琴图
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  # 只限制y轴（Naiveness分数）
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) + 
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  xlab('Naiveness score') +
  theme(legend.position = "none") +
  annotate("text", x = x_mid, y = y_max, label = paste0("p = ", signif(exhausted_pvalue, digits = 3)), size = 5, hjust = 0.5)+
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.1,
    color = "#2E8B57",
    size = 0.2
  )

p2 <- ggplot(scPDAC_T_CD4@meta.data, aes(x = group, y = treg, colour = group, fill = group)) +
  geom_violin(size = 1, alpha = 0.3) +  # 使用小提琴图
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  # 只限制y轴（Treg分数）
  scale_colour_manual(values = c("#5CACEE", "#E95D69")) + 
  scale_fill_manual(values = c("#5CACEE", "#E95D69")) +  
  xlab('Treg score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  annotate("text", x = x_mid, y = y_max, label = paste0("p = ", signif(exhausted_pvalue, digits = 3)), size = 5, hjust = 0.5)+
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.1,
    color = "#2E8B57",
    size = 0.2
  )

p <- p1 + p2
ggsave('./Result/Tcells/Naiveness_Treg_score_violin_CD4.pdf', plot = p, width = 7.86, height = 2.67, dpi = 300)


scPDAC_T_CD4@meta.data$group_label <- ifelse(scPDAC_T_CD4@meta.data$celltype == "CD4 Treg", "CD4 Treg", "Other")
treg_pvalue <- wilcox.test(treg ~ group_label, data = scPDAC_T_CD4@meta.data)$p.value
print(treg_pvalue)
scPDAC_T_CD4@meta.data$group_label <- ifelse(scPDAC_T_CD4@meta.data$celltype == "CD4 Naive", "CD4 Naive", "Other")
naive_pvalue <- wilcox.test(treg ~ group_label, data = scPDAC_T_CD4@meta.data)$p.value
print(naive_pvalue)

p1 <- ggplot(scPDAC_T_CD4@meta.data, aes(x = celltype, y = naiveness, colour = celltype, fill = celltype)) +
  geom_violin(size = 1, alpha = 0.3) + 
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14)) +
  xlab('Naiveness score') +
  theme(legend.position = "none")
p2 <- ggplot(scPDAC_T_CD4@meta.data, aes(x = celltype, y = treg, colour = celltype, fill = celltype)) +
  geom_violin(size = 1, alpha = 0.3) +  
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2, 2)) +  
  xlab('Treg score') +
  theme(strip.background = element_rect(color = "white", fill = "white", linewidth = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 14))
p <- p1 + p2
ggsave('./Result/Tcells/Naiveness_Treg_score_violin_CD4_celltype.tiff', plot = p, width = 12.86, height = 2.67, dpi = 300)


library(ggplot2)
library(reshape2)
library(scales)
high_var_genes_heatmap_data <- read.csv("./Data/Tcells/Gene_highly_variable_pseudotime.csv")
high_var_genes_heatmap_data <- as.data.frame(t(high_var_genes_heatmap_data))
dim(high_var_genes_heatmap_data)
head(high_var_genes_heatmap_data)

high_var_genes_heatmap_data_head <- head(high_var_genes_heatmap_data, 200)
high_var_genes_heatmap_data_head <- as.data.frame(t(apply(high_var_genes_heatmap_data_head, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
})))

high_var_genes_heatmap_data_top6000cols <- high_var_genes_heatmap_data_head[, 1:6000]
gene_mean_values <- rowMeans(high_var_genes_heatmap_data_top6000cols)
sorted_genes_by_mean <- names(sort(gene_mean_values, decreasing = TRUE))
top_60_genes <- sorted_genes_by_mean[1:75]
high_var_genes_heatmap_data_head <- high_var_genes_heatmap_data_head[top_60_genes, ]

high_var_genes_heatmap_data_last8000 <- high_var_genes_heatmap_data_head[, (ncol(high_var_genes_heatmap_data_head)-7999):ncol(high_var_genes_heatmap_data_head)]
gene_mean_values_last8000 <- rowMeans(high_var_genes_heatmap_data_last8000)
sorted_genes_by_mean_last8000 <- names(sort(gene_mean_values_last8000, decreasing = TRUE))
top_50_genes_last8000 <- sorted_genes_by_mean_last8000[26:75]
high_var_genes_heatmap_data_head <- high_var_genes_heatmap_data_head[top_50_genes_last8000, ]

high_var_genes_heatmap_data_head$Gene <- rownames(high_var_genes_heatmap_data_head)
high_var_genes_heatmap_data_head_melted <- melt(high_var_genes_heatmap_data_head, id.vars = "Gene")
colnames(high_var_genes_heatmap_data_head_melted) <- c("Gene", "Cell", "Expression")

library(zoo)
library(stats)

kernel <- rep(1/1000, 1000)

high_var_genes_heatmap_data_head_melted_2 <- high_var_genes_heatmap_data_head_melted %>%
  group_by(Gene) %>%
  mutate(Smoothed_Expression = stats::filter(Expression, kernel, sides = 2)) %>%
  mutate(Smoothed_Expression = na.locf(Smoothed_Expression, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(!is.na(Smoothed_Expression))  

p.heatmap_head_smoothed <- ggplot(high_var_genes_heatmap_data_head_melted_2, aes(x = Cell, y = Gene, fill = Smoothed_Expression)) +
  geom_tile() +
  scale_fill_gradientn(colors = c('#7895C1', '#A8CBDF', '#A8CBDF', '#F7FBC9', '#F5EBAE', '#F0C284', '#EF8B67')) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.title = element_blank(),  
        panel.grid = element_blank(),
        legend.position = "none") + 
  labs(title = NULL)  

ggsave('./Result/Tcells/genes_heatmap_pseudotime_smoothed.pdf', plot = p.heatmap_head_smoothed, width = 6, height = 4)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(DOSE)
library(biomaRt)

high_var_genes <- rownames(high_var_genes_heatmap_data_head) 

GOS <- NULL
eg <- bitr(high_var_genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
go <- enrichGO(eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='BP',
               pAdjustMethod = 'BH',
               pvalueCutoff = 1, 
               qvalueCutoff = 1,
               keyType = 'ENTREZID')
go <- setReadable(go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
deggo <- go@result
bgcount <- as.numeric(unlist(strsplit(go@result[1,3],'/'))[2])
deggo$GeneRatios <- deggo$Count / bgcount
GOS <- rbind(GOS, deggo)

GOS <- GOS[GOS$pvalue<=0.05 & GOS$GeneRatios>=0.01,]
# write.csv(GOS,'D:/Bio/PDAC/PDAC/Result/Tcells/GO_high_var_genes_upregulated.csv',quote = F)

top10_GO <- GOS[order(GOS$qvalue),][1:10,]
top10_GO$neg_log10_qvalue <- -log10(top10_GO$qvalue)
p.up <- ggplot(top10_GO, aes(x = reorder(Description, -neg_log10_qvalue), y = neg_log10_qvalue)) + 
  geom_bar(stat = "identity", fill = '#EF8B67', color = NA) + 
  theme_minimal() + 
  coord_flip() + 
  scale_y_reverse() +  
  labs(x = "GO Terms", 
       y = "-log10(q-value)") +  
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),  
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        plot.title = element_blank(), 
        panel.grid = element_blank(), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank())  

ggsave('./Result/Tcells/go_upregulated.pdf', plot = p.up, width = 6, height = 4)

scCD4 <- readRDS("./Data/Tcells/scCD4.rds") 
high_var_genes_markers_df <- FindMarkers(scCD4, group.by = "group", ident.1 = "Resist", ident.2 = "Sensitive", features = high_var_genes, logfc.threshold = 0.25)
high_var_genes_markers_df$SYMBOL <- rownames(high_var_genes_markers_df)
head(high_var_genes_markers_df)
#eg<- bitr(high_var_genes, fromType ="SYMBOL", toType =C("ENTREZID"), OrgDb ="org.Hs.eg.db")

ranked_genes <- high_var_genes_markers_df$avg_log2FC
names(ranked_genes) <- rownames(high_var_genes_markers_df)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

ranked_genes_entrez <- bitr(names(ranked_genes), 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
names(ranked_genes) <- ranked_genes_entrez$SYMBOL
gsea_result <- gseGO(geneList = ranked_genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL", 
                     ont = "BP",          
                     pvalueCutoff = 1,  
                     verbose = TRUE)


p <- ridgeplot(gsea_result,
               showCategory = 20,
               fill = "pvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
               core_enrichment = TRUE,#是否只使用 core_enriched gene
               label_format = 30,
               orderBy = "NES",
               decreasing = FALSE
)+
  theme(axis.text.y = element_text(size=8))
ggsave('./Result/Tcells/kegg_ridgeplot_upregulated.pdf',plot = p, width = 6, height = 8)



high_var_genes_heatmap_data_tail <- tail(high_var_genes_heatmap_data, 800)
high_var_genes_heatmap_data_tail <- as.data.frame(t(apply(high_var_genes_heatmap_data_tail, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
})))

high_var_genes_heatmap_data_tail_last5000 <- high_var_genes_heatmap_data_tail[, (ncol(high_var_genes_heatmap_data_tail)-4999):ncol(high_var_genes_heatmap_data_tail)]
gene_mean_values_tail_last5000 <- rowMeans(high_var_genes_heatmap_data_tail_last5000)
sorted_genes_by_mean_tail_last5000 <- names(sort(gene_mean_values_tail_last5000, decreasing = TRUE))
top_50_genes_tail_last5000 <- sorted_genes_by_mean_tail_last5000[1:150]
high_var_genes_heatmap_data_top50_tail_last5000 <- high_var_genes_heatmap_data_tail[top_50_genes_tail_last5000, ]
high_var_genes_heatmap_data_tail <- high_var_genes_heatmap_data_top50_tail_last5000

high_var_genes_heatmap_data_tail_top6000cols <- high_var_genes_heatmap_data_tail[, 1:6000]
gene_mean_values_tail_top6000cols <- rowMeans(high_var_genes_heatmap_data_tail_top6000cols)
sorted_genes_by_mean_tail_top6000cols <- names(sort(gene_mean_values_tail_top6000cols, decreasing = TRUE))
top_50_genes_tail_top6000cols <- sorted_genes_by_mean_tail_top6000cols[101:150]
high_var_genes_heatmap_data_top50_tail_top6000cols <- high_var_genes_heatmap_data_tail[top_50_genes_tail_top6000cols, ]
high_var_genes_heatmap_data_tail <- high_var_genes_heatmap_data_top50_tail_top6000cols

high_var_genes_heatmap_data_tail$Gene <- rownames(high_var_genes_heatmap_data_tail)
high_var_genes_heatmap_data_tail_melted <- melt(high_var_genes_heatmap_data_tail, id.vars = "Gene")
colnames(high_var_genes_heatmap_data_tail_melted) <- c("Gene", "Cell", "Expression")

kernel <- rep(1/1000, 1000)

high_var_genes_heatmap_data_tail_melted_2 <- high_var_genes_heatmap_data_tail_melted %>%
  group_by(Gene) %>%
  mutate(Smoothed_Expression = stats::filter(Expression, kernel, sides = 2)) %>%
  mutate(Smoothed_Expression = na.locf(Smoothed_Expression, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(!is.na(Smoothed_Expression))  

p.heatmap_upregrated <- ggplot(high_var_genes_heatmap_data_tail_melted_2, aes(x = Cell, y = Gene, fill = Smoothed_Expression)) +
  geom_tile() +
  scale_fill_gradientn(colors = c('#7895C1', '#A8CBDF', '#A8CBDF', '#F7FBC9', '#F5EBAE', '#F0C284', '#EF8B67')) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.title = element_blank(),  
        panel.grid = element_blank(),
        legend.position = "none") +  
  labs(title = NULL) 

ggsave('./Result/Tcells/genes_heatmap_pseudotime_upregrated.pdf', plot = p.heatmap_upregrated, width = 6, height = 4)



high_var_genes <- rownames(high_var_genes_heatmap_data_tail) 

GOS <- NULL
eg <- bitr(high_var_genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
go <- enrichGO(eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='BP',
               pAdjustMethod = 'BH',
               pvalueCutoff = 1, 
               qvalueCutoff = 1,
               keyType = 'ENTREZID')
go <- setReadable(go, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
deggo <- go@result
bgcount <- as.numeric(unlist(strsplit(go@result[1,3],'/'))[2])
deggo$GeneRatios <- deggo$Count / bgcount
GOS <- rbind(GOS, deggo)

GOS <- GOS[GOS$pvalue<=0.05 & GOS$GeneRatios>=0.01,]
write.csv(GOS,'D:/Bio/PDAC/PDAC/Result/Tcells/GO_high_var_genes_downregulated.csv',quote = F)

top10_GO <- GOS[order(GOS$qvalue),][1:10,]
top10_GO$neg_log10_qvalue <- -log10(top10_GO$qvalue)
p.down <- ggplot(top10_GO, aes(x = reorder(Description, -neg_log10_qvalue), y = neg_log10_qvalue)) + 
  geom_bar(stat = "identity", fill = '#A8CBDF', color = NA) +  # 去掉条形图边框颜色
  theme_minimal() + 
  coord_flip() +  
  scale_y_reverse() +  
  labs(x = "GO Terms", 
       y = "-log10(q-value)") + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_blank(), 
        panel.grid = element_blank(), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5),
        axis.ticks.y = element_blank())  

ggsave('./Result/Tcells/go_downregulated.pdf',plot = p.down, width = 6, height = 4)


scCD4 <- readRDS("./Data/Tcells/scCD4.rds") 
high_var_genes <- high_var_genes[high_var_genes %in% rownames(scCD4)]
high_var_genes_markers_df <- FindMarkers(scCD4, group.by = "group", ident.1 = "Resist", ident.2 = "Sensitive", features = high_var_genes, logfc.threshold = 0.25)
high_var_genes_markers_df$SYMBOL <- rownames(high_var_genes_markers_df)
head(high_var_genes_markers_df)
#eg<- bitr(high_var_genes, fromType ="SYMBOL", toType =C("ENTREZID"), OrgDb ="org.Hs.eg.db")

ranked_genes <- high_var_genes_markers_df$avg_log2FC
names(ranked_genes) <- rownames(high_var_genes_markers_df)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

ranked_genes_entrez <- bitr(names(ranked_genes), 
                            fromType = "SYMBOL", 
                            toType = "ENTREZID", 
                            OrgDb = org.Hs.eg.db)
names(ranked_genes) <- ranked_genes_entrez$SYMBOL

gsea_result <- gseGO(geneList = ranked_genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL", 
                     ont = "BP",          
                     pvalueCutoff = 1,  
                     verbose = TRUE)


p <- ridgeplot(gsea_result,
               showCategory = 20,
               fill = "pvalue", #填充色 "pvalue", "p.adjust", "qvalue" 
               core_enrichment = TRUE,#是否只使用 core_enriched gene
               label_format = 30,
               orderBy = "NES",
               decreasing = FALSE
)+
  theme(axis.text.y = element_text(size=8))

ggsave('./Result/Tcells/kegg_ridgeplot_downregulated.pdf',plot = p, width = 6, height = 8)


