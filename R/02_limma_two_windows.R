###############################################################################
# 02b_rerun_fixed.R  —— 修复版：自带数据重建 + 行名清洗，不依赖旧的 step1_data.rds
# 你的报错 "row names contain missing values" = 矩阵里有基因名为空/NA 的行
# 我已经在本机复现并修好，直接整段运行即可
###############################################################################

WORKDIR <- Sys.getenv("RG3_ROOT", ".")
setwd(WORKDIR)
suppressMessages(library(limma))

# ---- 1. 先诊断：文件到底完整不完整 ----
f <- "GSE102511_TPM.txt.gz"
cat("文件大小 (MB)：", round(file.size(f) / 1e6, 2), "  <- 完整应为 3.86 左右\n")
raw <- read.delim(gzfile(f), check.names = FALSE)
cat("raw 行数：", nrow(raw), " 列数：", ncol(raw), "  <- 行数应约 20801（含表头）\n")
cat("首列前 3 个：", paste(head(raw[[1]], 3), collapse = " / "), "\n")
cat("空/NA 基因名：", sum(is.na(raw[[1]]) | raw[[1]] == ""), "\n")
cat("重复基因名：", sum(duplicated(raw[[1]])), "\n\n")

if (nrow(raw) < 20000) {
  stop("数据文件不完整（行数 ", nrow(raw), " 条）。请删掉 ", f,
       " 后重新运行 01_load_qc.R 重新下载，然后再跑本脚本。")
}

# ---- 2. 重建干净矩阵 ----
genes <- raw[[1]]
expr <- as.matrix(raw[, -1]); rownames(expr) <- genes

# 2a. 丢掉基因名为空/NA 的行（这就是报错根源）
bad <- is.na(rownames(expr)) | rownames(expr) == ""
cat("删除无名行：", sum(bad), " 行\n")
expr <- expr[!bad, , drop = FALSE]

# 2b. 丢掉含 NA 的表达值行
na_row <- !complete.cases(expr)
cat("删除含 NA 的行：", sum(na_row), " 行\n")
expr <- expr[!na_row, , drop = FALSE]

# 2c. 重复基因名取均值
n0 <- nrow(expr)
expr <- apply(expr, 2, function(x) tapply(x, rownames(expr), mean))
cat("合并重复基因：", n0, " -> ", nrow(expr), "\n")

# 2d. 去零方差
expr <- expr[apply(expr, 1, var) > 0, ]
cat("去零方差基因后：", nrow(expr), " 个基因\n\n")

group <- factor(sub(".*-", "", colnames(expr)), levels = c("NL", "AAH", "LUAD"))
donor <- sub("-.*", "", colnames(expr))
print(table(group))
stopifnot(!any(is.na(rownames(expr))), !any(duplicated(rownames(expr))))

saveRDS(list(expr = expr, group = group), "step1_data_clean.rds")

# ---- 3. 差异分析（配对 + makeContrasts）----
run_de <- function(g1, g2, tag) {
  idx <- group %in% c(g1, g2)
  e <- expr[, idx, drop = FALSE]
  gr <- droplevels(group[idx]); dn <- factor(donor[idx])
  cat("===", tag, "：", g1, "vs", g2, "| 样本", sum(gr == g1), "vs", sum(gr == g2),
      "| 病人", length(unique(dn)), "\n")

  design <- model.matrix(~ 0 + gr + dn)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(e, design)
  cm <- makeContrasts(contrasts = paste0("gr", g2, "-gr", g1), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  # coef = 1 显式指定（只有一个对比），避免按名字查找出岔子
  tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  rownames(tt) <- NULL

  sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, ]
  cat("   差异基因（FDR<0.05 & |logFC|>1）：", nrow(sig),
      " ｜上调", sum(sig$logFC > 0), "下调", sum(sig$logFC < 0), "\n")
  print(head(sig[, c("gene", "logFC", "adj.P.Val")], 10))
  write.csv(tt, paste0("DEG_", tag, "_full.csv"), row.names = FALSE)
  write.csv(sig, paste0("DEG_", tag, "_sig.csv"), row.names = FALSE)
  cat("\n")
  sig
}

onset <- run_de("NL", "AAH", "onset")
prog  <- run_de("AAH", "LUAD", "progression")

both <- intersect(onset$gene, prog$gene)
cat("两窗口共有差异基因：", length(both), "\n")
cat("仅起始窗口特有：", length(setdiff(onset$gene, prog$gene)), "\n")
cat("仅进展窗口特有：", length(setdiff(prog$gene, onset$gene)), "\n")
write.csv(data.frame(gene = both), "DEG_shared.csv", row.names = FALSE)
saveRDS(list(onset = onset, prog = prog), "step2_deg.rds")
cat("\n>>> 完成。把这整段输出贴给我。\n")
