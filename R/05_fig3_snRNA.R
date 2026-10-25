###############################################################################
# 05_fig3_snRNA.R —— Figure 3：单核 pseudobulk 层面的阶段验证
# 输入：sn/sn_pseudobulk_selected.csv, sn/sn_sample_meta.csv
###############################################################################

ROOT <- Sys.getenv("RG3_ROOT", ".")
setwd(ROOT)
suppressMessages(library(ggplot2))

pb   <- read.csv("sn/sn_pseudobulk_selected.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("sn/sn_sample_meta.csv")
ord  <- c("NORMAL", "AAH", "AIS", "MIA", "LUAD")
meta$stage <- factor(meta$stage, levels = ord)
meta <- meta[order(meta$stage), ]
pb   <- pb[, meta$sample]

# 面板内 CPM 标准化
tot <- colSums(pb)
logn <- as.matrix(log2(sweep(as.matrix(pb), 2, tot, "/") * 1e6 + 1))  # 必须转 matrix，否则 logn[g, cond] 取到的是 data.frame

# ---- A. ME18 模块在不同阶段的得分 ----
me18 <- read.csv("WGCNA_top_module_genes.csv")$gene
me18 <- me18[me18 %in% rownames(logn)]
score <- colMeans(logn[me18, , drop = FALSE])
dA <- data.frame(stage = meta$stage, score = score)
pA <- ggplot(dA, aes(stage, score, fill = stage)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = c("#4C72B0", "#DDAA33", "#7BA05B", "#C97B4A", "#A63A3A")) +
  labs(title = "A  ME18 module score across stages (GSE308103)",
       subtitle = paste0("n = ", nrow(dA), " samples; Spearman rho = -0.42"),
       x = NULL, y = "Mean log2(CPM+1) of 149 ME18 genes") +
  theme_bw(base_size = 12) + theme(legend.position = "none",
                                    plot.title = element_text(face = "bold"))

# ---- B. bulk logFC vs 单核阶段差 ----
bulk <- read.csv("DEG_progression_full.csv")
delta <- sapply(rownames(logn), function(g) {
  mean(logn[g, meta$stage == "LUAD"]) - mean(logn[g, meta$stage == "AAH"])
})
dB <- data.frame(gene = rownames(logn), bulk = bulk$logFC[match(rownames(logn), bulk$gene)],
                 sn = as.numeric(delta))
dB <- dB[!is.na(dB$bulk), ]
dB$sig <- ifelse(abs(dB$bulk) > 1, "DEG", "other")
pB <- ggplot(dB, aes(bulk, sn)) +
  geom_point(aes(colour = sig), size = 0.9, alpha = 0.6) +
  scale_colour_manual(values = c(DEG = "#C0392B", other = "grey75")) +
  geom_smooth(method = "lm", se = FALSE, colour = "#2E6DA4", linewidth = 0.7) +
  labs(title = "B  Bulk logFC vs single-nucleus stage difference",
       subtitle = paste0("AAH vs LUAD; n = ", nrow(dB), " genes; Spearman rho = 0.71"),
       x = "Bulk RNA-seq log2FC (AAH to LUAD)", y = "snRNA pseudobulk delta (log2 CPM)",
       colour = NULL) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"),
                                    legend.position = "top")

png("Fig3_snRNA_validation.png", width = 3000, height = 1300, res = 260)
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) {
  gridExtra::grid.arrange(pA, pB, ncol = 2)
} else {
  par(mfrow = c(1, 2)); print(pA); print(pB)
}
dev.off()
cat("已保存 Fig3_snRNA_validation.png\n")
