rm(list = ls())

library(survival)
library(survminer)
library(GEOquery)
library(pROC)
library(ggplot2)
library(gridExtra)

gse <- getGEO("GSE28735", GSEMatrix = TRUE)
#gse <- getGEO(filename = "D:/Bio/PDAC/PDAC/Code/Data/survival/GSE28735_series_matrix.txt.gz", GSEMatrix = TRUE)


expr_matrix <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])
feature_data <- fData(gse[[1]])

survival_time <- as.numeric(pheno_data$`survival_month:ch1`) * 30.44
survival_status <- as.numeric(pheno_data$`cancer_death:ch1`)
tissue_type <- pheno_data$`tissue:ch1`

meta <- data.frame(
  SampleID = rownames(pheno_data),
  time = survival_time,
  event = survival_status,
  tissue = tissue_type
)

meta <- meta[meta$tissue == "T", ]
meta <- na.omit(meta)

genes <- c("SPP1", "CD44")
gene_symbols <- feature_data$gene_assignment
probe_ids <- c()

for(gene in genes) {
  probe_id <- rownames(feature_data)[grep(gene, gene_symbols, ignore.case = TRUE)]
  if(length(probe_id) > 0) {
    probe_ids <- c(probe_ids, probe_id[1])
    print(paste("找到基因", gene, "的探针:", probe_id[1]))
  } else {
    print(paste("未找到基因", gene, "的探针"))
  }
}

if(length(probe_ids) > 0) {
  expr_gene <- expr_matrix[probe_ids, , drop = FALSE]
  
  common_samples <- intersect(colnames(expr_gene), meta$SampleID)
  expr_gene <- expr_gene[, common_samples, drop = FALSE]
  meta <- meta[meta$SampleID %in% common_samples, ]
  
  expr_gene <- t(expr_gene)
  colnames(expr_gene) <- genes
  
  dat <- cbind(meta, expr_gene)
  print(paste("最终分析样本数:", nrow(dat)))
  
  dir.create("./Result/fig7/GSE28735", showWarnings = FALSE, recursive = TRUE)
  
  dat$SPP1_CD44 <- dat$SPP1 * dat$CD44
  
  font_size <- 14
  axis_text_size <- 18
  axis_title_size <- 18
  pval_size <- 18
  
  large_theme <- theme(
    text = element_text(size = font_size),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    axis.title.y = element_text(size = axis_title_size),
    plot.title = element_text(size = font_size + 2, face = "bold"),
    legend.text = element_text(size = font_size),
    legend.title = element_text(size = font_size),
    strip.text = element_text(size = font_size)
  )
  
  increase_pval_size <- function(plot_obj, size) {
    for (i in 1:length(plot_obj$layers)) {
      layer <- plot_obj$layers[[i]]
      if ("GeomText" %in% class(layer$geom) || "GeomLabel" %in% class(layer$geom)) {
        if (!is.null(layer$aes_params)) {
          layer$aes_params$size <- size
        } else {
          layer$aes_params <- list(size = size)
        }
        if (!is.null(layer$mapping) && "size" %in% names(layer$mapping)) {
          layer$mapping$size <- NULL
        }
        plot_obj$layers[[i]] <- layer
      }
    }
    return(plot_obj)
  }
  
    roc_spp1_optimal <- roc(dat$event, dat$SPP1, direction = "<")
    optimal_cutoff_spp1 <- coords(roc_spp1_optimal, "best", ret = "threshold", best.method = "youden")$threshold
    cat("SPP1 最优截断点:", optimal_cutoff_spp1, "\n")
    
    dat$SPP1_optimal_group <- ifelse(dat$SPP1 > optimal_cutoff_spp1, "Low", "High")
    fit_spp1_optimal <- survfit(Surv(time, event) ~ SPP1_optimal_group, data = dat)
    p_spp1_optimal <- ggsurvplot(fit_spp1_optimal, data = dat,
                   palette=c("#00AFBB", "#FC4E07"),
                 pval=TRUE, risk.table=TRUE,
                   risk.table.height=0.25, risk.table.fontsize=4,
                 title="SPP1 expression and survival in GSE28735",
                   legend.labs=c("Low", "High"))
  p_spp1_optimal$plot <- p_spp1_optimal$plot + large_theme
  p_spp1_optimal$table <- p_spp1_optimal$table + large_theme + theme(legend.position = "none")
  p_spp1_optimal$plot <- increase_pval_size(p_spp1_optimal$plot, pval_size / .pt)
  combined_spp1 <- grid.arrange(p_spp1_optimal$plot, p_spp1_optimal$table, nrow = 2, heights = c(0.75, 0.25))
  ggsave("./Result/fig7/GSE28735/GSE28735_SPP1_optimal_cutoff_survival.pdf", plot = combined_spp1, width = 6, height = 8, device = "pdf")
  
    roc_cd44_optimal <- roc(dat$event, dat$CD44, direction = "<")
    optimal_cutoff_cd44 <- coords(roc_cd44_optimal, "best", ret = "threshold", best.method = "youden")$threshold
    cat("CD44 最优截断点:", optimal_cutoff_cd44, "\n")
    
    dat$CD44_optimal_group <- ifelse(dat$CD44 > optimal_cutoff_cd44, "Low", "High")
    fit_cd44_optimal <- survfit(Surv(time, event) ~ CD44_optimal_group, data = dat)
    p_cd44_optimal <- ggsurvplot(fit_cd44_optimal, data = dat,
                   palette=c("#00AFBB", "#FC4E07"),
                 pval=TRUE, risk.table=TRUE,
                   risk.table.height=0.25, risk.table.fontsize=4,
                 title="CD44 expression and survival in GSE28735",
                   legend.labs=c("Low", "High"))
  p_cd44_optimal$plot <- p_cd44_optimal$plot + large_theme
  p_cd44_optimal$table <- p_cd44_optimal$table + large_theme + theme(legend.position = "none")
  p_cd44_optimal$plot <- increase_pval_size(p_cd44_optimal$plot, pval_size / .pt)
  combined_cd44 <- grid.arrange(p_cd44_optimal$plot, p_cd44_optimal$table, nrow = 2, heights = c(0.75, 0.25))
  ggsave("./Result/fig7/GSE28735/GSE28735_CD44_optimal_cutoff_survival.pdf", plot = combined_cd44, width = 6, height = 8, device = "pdf")
  
    roc_combined_optimal <- roc(dat$event, dat$SPP1_CD44, direction = "<")
    optimal_cutoff_combined <- coords(roc_combined_optimal, "best", ret = "threshold", best.method = "youden")$threshold
    cat("SPP1×CD44 最优截断点:", optimal_cutoff_combined, "\n")
    
    dat$SPP1_CD44_optimal_group <- ifelse(dat$SPP1_CD44 > optimal_cutoff_combined, "Low", "High")
    fit_combined_optimal <- survfit(Surv(time, event) ~ SPP1_CD44_optimal_group, data = dat)
    p_combined_optimal <- ggsurvplot(fit_combined_optimal, data = dat,
                   palette=c("#00AFBB", "#FC4E07"),
                 pval=TRUE, risk.table=TRUE,
                   risk.table.height=0.25, risk.table.fontsize=4,
                 title="SPP1×CD44 signature in GSE28735",
                   legend.labs=c("Low", "High"))
  p_combined_optimal$plot <- p_combined_optimal$plot + large_theme
  p_combined_optimal$table <- p_combined_optimal$table + large_theme + theme(legend.position = "none")
  p_combined_optimal$plot <- increase_pval_size(p_combined_optimal$plot, pval_size / .pt)
  combined_optimal <- grid.arrange(p_combined_optimal$plot, p_combined_optimal$table, nrow = 2, heights = c(0.75, 0.25))
  ggsave("./Result/fig7/GSE28735/GSE28735_SPP1_CD44_optimal_cutoff_survival.pdf", plot = combined_optimal, width = 6, height = 8, device = "pdf")
  
    optimal_cutoffs <- data.frame(
      Marker = c("SPP1", "CD44", "SPP1×CD44"),
      Optimal_Cutoff = c(optimal_cutoff_spp1, optimal_cutoff_cd44, optimal_cutoff_combined),
      AUC = c(roc_spp1_optimal$auc, roc_cd44_optimal$auc, roc_combined_optimal$auc)
    )
    
    cat("\n最优截断点结果:\n")
    print(optimal_cutoffs)
  write.csv(optimal_cutoffs, "./Result/fig7/GSE28735/optimal_cutoffs.csv", row.names = FALSE)
  
} else {
  print("未找到目标基因，请检查基因名称或探针注释")
}
