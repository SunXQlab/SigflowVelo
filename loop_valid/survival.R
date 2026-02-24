library(survival)
library(survminer)
library(pROC)
library(ggplot2)
library(gridExtra)

proj='./loop_valid/tcga-paad'
load(file = paste0(proj,".for_survival.rdata") )

genes <- c("SPP1","CD44")
expr_gene <- t(exprSet[genes, rownames(meta)])  
colnames(expr_gene) <- genes

dat <- cbind(meta, expr_gene)
head(dat)

safe_pdf_plot <- function(filename, ggsurv_obj, width = 6, height = 8) {
  tryCatch({
    output_dir <- dirname(filename)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    if ("table" %in% names(ggsurv_obj)) {
      combined_plot <- grid.arrange(ggsurv_obj$plot, ggsurv_obj$table, 
                                   nrow = 2, heights = c(0.75, 0.25))
      ggsave(filename, plot = combined_plot, width = width, height = height, device = "pdf")
    } else {
      ggsave(filename, plot = ggsurv_obj$plot, width = width, height = height, device = "pdf")
    }
    return(TRUE)
  }, error = function(e) {
    tryCatch({
      pdf(filename, width = width, height = height, onefile = FALSE)
      print(ggsurv_obj)
      dev.off()
      return(TRUE)
    }, error = function(e2) {
      if (dev.cur() != 1) try(dev.off(), silent = TRUE)
      cat("保存失败:", filename, "错误:", e$message, "\n")
    return(FALSE)
    })
  })
}

dat$SPP1_CD44 <- dat$SPP1 * dat$CD44

roc_spp1_optimal <- roc(dat$event, dat$SPP1, direction = "<")
optimal_cutoff_spp1 <- coords(roc_spp1_optimal, "best", ret = "threshold", best.method = "youden")$threshold
cat("SPP1 最优截断点:", optimal_cutoff_spp1, "\n")

roc_cd44_optimal <- roc(dat$event, dat$CD44, direction = "<")
optimal_cutoff_cd44 <- coords(roc_cd44_optimal, "best", ret = "threshold", best.method = "youden")$threshold
cat("CD44 最优截断点:", optimal_cutoff_cd44, "\n")

roc_combined_optimal <- roc(dat$event, dat$SPP1_CD44, direction = "<")
optimal_cutoff_combined <- coords(roc_combined_optimal, "best", ret = "threshold", best.method = "youden")$threshold
cat("SPP1×CD44 最优截断点:", optimal_cutoff_combined, "\n")

dat$SPP1_optimal_group <- ifelse(dat$SPP1 > optimal_cutoff_spp1, "Low", "High")
dat$CD44_optimal_group <- ifelse(dat$CD44 > optimal_cutoff_cd44, "Low", "High")
dat$SPP1_CD44_optimal_group <- ifelse(dat$SPP1_CD44 > optimal_cutoff_combined, "Low", "High")

pval_size <- 16
axis_text_size <- 16

fit_spp1_optimal <- survfit(Surv(time, event) ~ SPP1_optimal_group, data = dat)
p_spp1_optimal <- ggsurvplot(fit_spp1_optimal, data = dat,
           palette=c("#00AFBB", "#FC4E07"),
           pval=TRUE, risk.table=TRUE,
           title="SPP1 expression and survival in PAAD",
           legend.labs=c("Low", "High"),
           fontsize = 5)

fit_cd44_optimal <- survfit(Surv(time, event) ~ CD44_optimal_group, data = dat)
p_cd44_optimal <- ggsurvplot(fit_cd44_optimal, data = dat,
           palette=c("#00AFBB", "#FC4E07"),
           pval=TRUE, risk.table=TRUE,
           title="CD44 expression and survival in PAAD",
           legend.labs=c("Low", "High"),
           fontsize = 5)

fit_combined_optimal <- survfit(Surv(time, event) ~ SPP1_CD44_optimal_group, data = dat)
p_combined_optimal <- ggsurvplot(fit_combined_optimal, data = dat,
           palette=c("#00AFBB", "#FC4E07"),
           pval=TRUE, risk.table=TRUE,
           title="SPP1×CD44 expression and survival in PAAD ",
           legend.labs=c("Low", "High"),
           fontsize = 5)

output_dir <- "./Result/fig7"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

combined_plots <- grid.arrange(
  p_spp1_optimal$plot, 
  p_cd44_optimal$plot, 
  p_combined_optimal$plot,
  p_spp1_optimal$table, 
  p_cd44_optimal$table, 
  p_combined_optimal$table,
  nrow = 2, ncol = 3,
  heights = c(0.75, 0.25)
)

ggsave("./Result/fig7/combined_survival_analysis.pdf", 
       plot = combined_plots, 
       width = 18, height = 8, device = "pdf")
cat("组合图已保存到: ./Result/fig7/combined_survival_analysis.pdf\n")

optimal_cutoffs <- data.frame(
  Marker = c("SPP1", "CD44", "SPP1×CD44"),
  Optimal_Cutoff = c(optimal_cutoff_spp1, optimal_cutoff_cd44, optimal_cutoff_combined),
  AUC = c(roc_spp1_optimal$auc, roc_cd44_optimal$auc, roc_combined_optimal$auc)
)

cat("最优截断点结果:\n")
print(optimal_cutoffs)
write.csv(optimal_cutoffs, "./Result/fig7/optimal_cutoffs.csv", row.names = FALSE)
