###############################################################################
# 01_load_qc.R  —— 第一枪 Step 1：载入 GSE102511 表达矩阵 + 质控
# 用法：把下面 WORKDIR 改成你自己的文件夹，然后整段选中运行（Ctrl+Enter）
###############################################################################

# ---- 0. 设定工作目录（改成你自己的；路径里的斜杠要用 / 或 \\）----
WORKDIR <- Sys.getenv("RG3_ROOT", ".")
if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE)
setwd(WORKDIR)
cat("工作目录：", getwd(), "\n")

# ---- 1. 下载数据（只需第一次；已有就跳过）----
f <- "GSE102511_TPM.txt.gz"
if (!file.exists(f)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102511/suppl/GSE102511_Smruthy_etal_allsamples_TPM.txt.gz"
  cat("正在下载（约 4 MB，耐心等）...\n")
  tryCatch(
    download.file(url, f, mode = "wb", method = "libcurl"),
    error = function(e) message("自动下载失败：", e$message,
      "\n请手动用浏览器打开上面的 url，把文件存到 ", WORKDIR)
  )
}
stopifnot(file.exists(f))

# ---- 2. 读取表达矩阵 ----
# 注意：文件里有重复基因名，不能直接用 row.names=1，否则报错
raw <- read.delim(gzfile(f), check.names = FALSE)
genes <- raw[[1]]
expr <- as.matrix(raw[, -1])
rownames(expr) <- genes
cat("矩阵维度：", nrow(expr), "个基因 x", ncol(expr), "个样本\n")
cat("数值范围：", round(min(expr), 2), "~", round(max(expr), 2),
    "（log2(TPM+1)，无需再取对数）\n")

# 重复基因名取平均值（必须在分组解析之前做）
if (any(duplicated(rownames(expr)))) {
  n0 <- nrow(expr)
  expr <- apply(expr, 2, function(x) tapply(x, rownames(expr), mean))
  cat("合并重复基因：", n0, "->", nrow(expr), "\n")
}

# ---- 3. 从列名解析分组 ----
# 列名形如 001-NL / 002-AAH / 003-LUAD
group <- sub(".*-", "", colnames(expr))
group <- factor(group, levels = c("NL", "AAH", "LUAD"))
print(table(group))
stopifnot(!any(is.na(group)))   # 有 NA 说明列名格式变了，立刻停下来看

# ---- 4. 清理 ----
# 4a. 去掉全零/全同值的基因
keep <- apply(expr, 1, function(x) sum(x > 0) >= 3)
expr <- expr[keep, ]
cat("减去低表达基因后：", nrow(expr), "个基因\n")

# 4b. 重复基因名取平均值（TPM 文件里经常有重名）
if (any(duplicated(rownames(expr)))) {
  n0 <- nrow(expr)
  expr <- apply(expr, 2, function(x) tapply(x, rownames(expr), mean))
  cat("合并重复基因：", n0, "->", nrow(expr), "\n")
}

# ---- 5. 保存 ----
saveRDS(list(expr = expr, group = group), "step1_data.rds")
write.csv(data.frame(sample = colnames(expr), group = as.character(group)),
          "sample_group.csv", row.names = FALSE)
cat("已保存 step1_data.rds 和 sample_group.csv\n")

# ---- 6. 质控图：PCA（看三组能不能分开）----
v <- apply(expr, 1, var)
top <- names(sort(v, decreasing = TRUE))[1:1000]
pca <- prcomp(t(expr[top, ]), scale. = TRUE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cols <- c("#2E7D46", "#C89B2C", "#8A3A3A")[as.numeric(group)]

png("QC_PCA.png", width = 1000, height = 800, res = 130)
plot(pca$x[, 1], pca$x[, 2], pch = 19, col = cols, cex = 1.4,
     xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
     main = "GSE102511: NL vs AAH vs LUAD")
legend("topright", legend = levels(group), col = c("#2E7D46", "#C89B2C", "#8A3A3A"),
       pch = 19, bty = "n")
dev.off()
cat("已保存 QC_PCA.png\n")
print(summary(pca$x[, 1:2]))
cat("\n>>> Step 1 完成。把屏幕上的全部输出（含矩阵维度、table(group)）贴给我。\n")
