###############################################################################
# 03_wgcna.R  —— Step 3A：基因名修复 + WGCNA（找与"结节阶段"相关的基因模块）
# 前提：已跑完 02b_rerun_fixed.R（需要 step1_data_clean.rds）
# 依赖：WGCNA（首次：BiocManager::install("WGCNA")，装得慢，等 10–20 分钟）
###############################################################################

WORKDIR <- Sys.getenv("RG3_ROOT", ".")
setwd(WORKDIR)
suppressMessages({library(WGCNA)})
options(stringsAsFactors = FALSE)
# macOS + RStudio 下开多线程常导致 "scheduled cores did not deliver results" 直接报错
# 这里强制单线程（慢一点但稳），遇到同样的错就靠这一行
disableWGCNAThreads()

# 读入数据：优先用 02b 生成的 step1_data_clean.rds；没有就自动从原始 TPM 重建
# （不管你之前跑的是 02 还是 02b，都不会再报"文件不存在"）
if (file.exists("step1_data_clean.rds")) {
  d <- readRDS("step1_data_clean.rds")
  expr <- d$expr; group <- d$group
  cat("载入已有 step1_data_clean.rds：", nrow(expr), "基因 x", ncol(expr), "样本\n")
} else {
  cat("未找到 step1_data_clean.rds，正在从 GSE102511_TPM.txt.gz 重建...\n")
  raw <- read.delim(gzfile("GSE102511_TPM.txt.gz"), check.names = FALSE)
  genes <- raw[[1]]
  expr <- as.matrix(raw[, -1]); rownames(expr) <- genes
  bad <- is.na(rownames(expr)) | rownames(expr) == ""
  expr <- expr[!bad, , drop = FALSE]
  expr <- expr[complete.cases(expr), , drop = FALSE]
  if (any(duplicated(rownames(expr))))
    expr <- apply(expr, 2, function(x) tapply(x, rownames(expr), mean))
  expr <- expr[apply(expr, 1, var) > 0, , drop = FALSE]
  group <- factor(sub(".*-", "", colnames(expr)), levels = c("NL", "AAH", "LUAD"))
  stopifnot(nrow(expr) > 15000)
  saveRDS(list(expr = expr, group = group), "step1_data_clean.rds")
  cat("重建完成：", nrow(expr), "基因 x", ncol(expr), "样本\n")
}

# ---- 1. 修复被 Excel 破坏的基因名（MARCH1 变成 1-Mar，SEPT12 变成 12-Sep）----
mon <- c(Jan = "JAN", Feb = "FEB", Mar = "MARCH", Apr = "APR", May = "MAY",
         Jun = "JUN", Jul = "JUL", Aug = "AUG", Sep = "SEPT", Oct = "OCT",
         Nov = "NOV", Dec = "DEC")
fix <- function(x) {
  i <- grepl("^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$", x)
  if (any(i)) x[i] <- paste0(mon[sub(".*-", "", x[i])], sub("-.*", "", x[i]))
  x
}
old <- rownames(expr)
new <- fix(old)
cat("修复的基因名：", sum(old != new), " 个\n")
print(head(data.frame(old = old[old != new], new = new[old != new]), 10))
rownames(expr) <- new

# ---- 2. 转置（WGCNA 要 样本 x 基因）----
dat <- t(expr)
gsg <- goodSamplesGenes(dat, verbose = 0)
cat("坏样本/坏基因：", sum(!gsg$goodSamples), "/", sum(!gsg$goodGenes), "\n")
dat <- dat[gsg$goodSamples, gsg$goodGenes]

# ---- 3. 选软阈值 ----
# 软阈值只在"变异最大的 5000 个基因"上估计（全基因算会吃掉几 GB 内存）
# 常规做法，Methods 里写明即可
powers <- c(1:10, seq(12, 20, 2))
v <- apply(dat, 2, var)
dat_soft <- dat[, names(sort(v, decreasing = TRUE))[1:min(5000, ncol(dat))]]
sft <- tryCatch(
  pickSoftThreshold(dat_soft, powerVector = powers, networkType = "signed", verbose = 0),
  error = function(e) { cat("pickSoftThreshold 失败：", conditionMessage(e),
                             "\n改用 power = 7（本数据已验证值）\n"); NULL }
)
tab <- if (is.null(sft)) 7 else sft$powerEstimate
if (!is.null(sft)) print(data.frame(power = sft$fitIndices$Power,
                 R2 = round(-sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq, 3)))
if (is.na(tab)) tab <- 7   # 没有达到 0.8 时退而用 7（本数据已验证，稿子里要说明）
cat("采用的软阈值 power =", tab, "\n")

# ---- 4. 建网络、找模块 ----
net <- blockwiseModules(dat, power = tab, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        maxBlockSize = 5000,
                        pamRespectsDendro = FALSE, verbose = 0)
modules <- net$colors
cat("共识别模块数：", length(unique(modules)) - (0 %in% modules), "\n")
print(table(modules))

# ---- 5. 模块与"疾病阶段"的关系（核心结果）----
stage <- as.numeric(group) - 1        # NL=0, AAH=1, LUAD=2
donor <- as.numeric(factor(sub("-.*", "", rownames(dat))))
MEs <- moduleEigengenes(dat, modules)$eigengenes
MEs <- orderMEs(MEs)
cor_stage <- cor(MEs, stage, use = "p")
p_stage   <- corPvalueStudent(cor_stage, nrow(dat))
cor_donor <- cor(MEs, donor, use = "p")

res <- data.frame(module = rownames(cor_stage),
                  r_stage = round(cor_stage[, 1], 3),
                  p_stage = signif(p_stage[, 1], 3),
                  r_donor = round(cor_donor[, 1], 3))
# 去掉 grey 模块（未分配基因）
res <- res[res$module != "MEgrey", ]
res <- res[order(-abs(res$r_stage)), ]
cat("\n=== 模块 vs 疾病阶段（按相关性排序）===\n")
print(res)
write.csv(res, "WGCNA_module_stage.csv", row.names = FALSE)

# 最相关模块的基因
top_mod <- res$module[1]
mod_num <- as.numeric(sub("ME", "", top_mod))
top_genes <- colnames(dat)[modules == mod_num]   # dat 是 样本x基因，基因名在列名上
cat("\n与阶段最相关的模块：", top_mod, "｜基因数：", length(top_genes), "\n")
write.csv(data.frame(gene = top_genes), "WGCNA_top_module_genes.csv", row.names = FALSE)

# 全部模块基因（供后续与 Rg3 靶点取交集）
allmap <- data.frame(gene = colnames(dat), module = paste0("ME", modules))
write.csv(allmap, "WGCNA_all_modules.csv", row.names = FALSE)

saveRDS(list(modules = modules, MEs = MEs, power = tab), "step3_wgcna.rds")

# ---- 6. 画一张模块-性状热图 ----
png("WGCNA_module_stage_heatmap.png", width = 900, height = 300 + 22 * nrow(res), res = 130)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = as.matrix(res$r_stage), xLabels = "Disease stage",
               yLabels = res$module, ySymbols = res$module,
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = res$r_stage, setStdMargins = FALSE,
               cex.text = 0.7, zlim = c(-1, 1),
               main = "Module-stage correlation (NL=0, AAH=1, LUAD=2)")
dev.off()
cat("\n已保存：WGCNA_module_stage.csv / WGCNA_top_module_genes.csv / WGCNA_all_modules.csv / WGCNA_module_stage_heatmap.png\n")
cat(">>> 把 '模块 vs 疾病阶段' 那张表贴给我。\n")
