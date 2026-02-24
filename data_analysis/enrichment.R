library(msigdbr)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)

scPDAC <- readRDS("./Data/scPDAC_anno.rds")

deglist <- list()
for(c in unique(scPDAC$celltype)){
  sce1 <- subset(scPDAC, celltype==c)
  degs <- FindMarkers(sce1,
                      ident.1 = 'Sensitive',
                      ident.2 = 'Resist',
                      assay = "RNA",
                      min.pct = 0.05,
                      logfc.threshold = 0.2,
                      group.by = 'group')
  degs <- cbind(gene=rownames(degs),degs)
  degs$celltype <- c
  deglist[[c]] <- degs
}
degs <- Reduce(rbind,deglist)
write.csv(degs,'./Data/celltype_group_degs.csv',quote = F,row.names = F)

degs <- read.csv('./Data/celltype_group_degs.csv')
setwd('./Result/Fig1/')

Cell_type <- unique(degs$celltype)
Cell_type_GO <- matrix(NA,0,0)

for (name in Cell_type) {
  dir.create(name, showWarnings = FALSE)
}

for (k in 1:length(Cell_type)){
  
  DEG.gene_symbol <- as.vector(as.character(as.matrix(degs$gene[which(degs$celltype==Cell_type[k])]))) 
  
  DEG.entrez_id <- mapIds(x = org.Hs.eg.db,
                         keys = DEG.gene_symbol,
                         keytype = "SYMBOL",
                         column = "ENTREZID")
  DEG.entrez_id <- na.omit(DEG.entrez_id)
  DEG.entrez_id <- data.frame(DEG.entrez_id)
  
  GO_enrich <- clusterProfiler::enrichGO(gene = DEG.entrez_id[,1],
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'ENTREZID',
                                         ont = 'ALL',
                                         pAdjustMethod = 'fdr',
                                         pvalueCutoff = 1,
                                         qvalueCutoff = 1,
                                         readable = FALSE)
  GO_enrich <- data.frame(GO_enrich)
  
  write.csv(GO_enrich, 
            paste(paste('./', Cell_type[k], sep = ''), 'GO_enrich.csv', sep = "/"),
            row.names = TRUE)
  
  Cell_type_GO <- rbind(Cell_type_GO, cbind(Cell_type[k], GO_enrich))
  
}

celltypes <- list.dirs(path = './',full.names = T,recursive = TRUE)[-1]
golist <- list()
for(c in celltypes){
  GO <- read.csv(paste(c,'/GO_enrich.csv',sep = ''),row.names = 1)
  GO$celltype <- gsub('./','',c)
  golist[[c]] <- GO
}
df1 <- Reduce(rbind,golist)

df1 <- df1[df1$p.adjust<=0.05,]
gobp <- df1[df1$ONTOLOGY=='BP',]
gobp$GeneRatio <- as.numeric(sub("(.*)/(.*)", "\\1", gobp$GeneRatio)) / as.numeric(sub("(.*)/(.*)", "\\2", gobp$GeneRatio))
gobp <- gobp %>% group_by(celltype) %>% arrange(desc(GeneRatio)) %>% do(head(., n = 10))

paths = c('regulation of T cell activation', 
           'immune response-activating signaling pathway',
           'immune response-regulating signaling pathway',
           'leukocyte mediated immunity',
           'leukocyte cell-cell adhesion',
           'regulation of inflammatory response',
           'leukocyte migration',
           'myeloid cell differentiation',
           'chemotaxis',
           'epithelial cell proliferation',
           'cell growth',
           'small GTPase mediated signal transduction',
           'positive regulation of cytokine production')
gobp <- filter(gobp, Description %in% paths)
gobp$Description <- factor(gobp$Description,levels = unique(gobp$Description))
gobp$celltype <- factor(gobp$celltype, levels = c("Acinar", "B cells", "T cells", "NK cells", "Dutcal", "Endothelial", "Fibroblast", "Macrophage", "Mast", "Plasma", "Neutroohils"))
p.gobp.generatio <- ggplot(gobp,aes(x=celltype,y=Description)) +
  geom_point(aes(size=GeneRatio,fill=p.adjust),shape = 21,stroke = 1) +
  scale_fill_gradient(high='#336699',low = "#FF9896FF") +
  theme_bw(base_size = 18) +
  theme(text=element_text(size=16,  family="sans"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(0.01, "cm")) +
  ggtitle('GO-BP') +
  xlab(NULL) +
  ylab(NULL)

ggsave('celltype_GOBP_generatio_selected.pdf',plot = p.gobp.generatio,width = 9,height = 5)
