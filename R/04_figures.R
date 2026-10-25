###############################################################################
# 04_figures.R  —— Step 4：出图（火山图 ×2、桥接基因热图、PPI 网络图）
# 前提：已跑完 02（或 02b）、03_wgcna.R，并已放好 rg3_targets_swisstarget.csv
# 依赖：ggplot2、igraph
# 说明：桥接基因自己算；PPI 边缺失时自动从 STRING 下载（不需要我这边传文件）
###############################################################################

WORKDIR <- Sys.getenv("RG3_ROOT", ".")
setwd(WORKDIR)
suppressMessages({library(ggplot2); library(igraph)})

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold", size = 13))

# ============ 0. 读靶点 + 自己算桥接基因 ============
if (!file.exists("rg3_targets_swisstarget.csv"))
  stop("找不到 rg3_targets_swisstarget.csv，请把靶点表放到：", normalizePath(WORKDIR))
tg <- read.csv("rg3_targets_swisstarget.csv")
tg <- unique(tg$gene[!tg$gene %in% c("N/A", "N", "")])
cat("Rg3 靶点数：", length(tg), "\n")

onset_f <- read.csv("DEG_onset_full.csv")
prog_f  <- read.csv("DEG_progression_full.csv")
sig_onset <- onset_f$gene[onset_f$adj.P.Val < 0.05 & abs(onset_f$logFC) > 0.5]
sig_prog  <- prog_f$gene[prog_f$adj.P.Val  < 0.05 & abs(prog_f$logFC)  > 0.5]
bridge <- sort(intersect(tg, c(sig_onset, sig_prog)))
cat("桥接基因（靶点 ∩ 显著差异基因）：", length(bridge), "个\n")
write.csv(data.frame(gene = bridge), "bridge_genes.csv", row.names = FALSE)

# ============ 1. 火山图 ============
volcano <- function(file, title, out, lfc = 1) {
  d <- read.csv(file)
  d$sig <- ifelse(d$adj.P.Val < 0.05 & abs(d$logFC) > lfc,
                  ifelse(d$logFC > 0, "Up", "Down"), "NS")
  d$sig <- factor(d$sig, levels = c("Up", "Down", "NS"))
  cat(title, "：上调", sum(d$sig == "Up"), "下调", sum(d$sig == "Down"), "\n")

  lab <- rbind(head(d[d$sig == "Up", ][order(-d[d$sig == "Up", ]$logFC), ], 6),
               head(d[d$sig == "Down", ][order(d[d$sig == "Down", ]$logFC), ], 6))
  lab <- lab[!duplicated(lab$gene), ]

  p <- ggplot(d, aes(logFC, -log10(adj.P.Val), colour = sig)) +
    geom_point(size = 1.1, alpha = 0.75) +
    scale_colour_manual(values = c(Up = "#C0392B", Down = "#2E6DA4", NS = "grey78")) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
    geom_text(data = lab, aes(label = gene), size = 3, vjust = -0.6, show.legend = FALSE) +
    labs(title = title, x = expression(log[2]~"fold change"), y = expression(-log[10]~"FDR"),
         colour = NULL) +
    theme_pub + theme(legend.position = "top")
  ggsave(out, p, width = 6.2, height = 5.4, dpi = 300)
  cat("  已保存", out, "\n")
}

volcano("DEG_onset_full.csv", "Normal lung vs AAH (onset window)", "Fig2A_volcano_onset.png")
volcano("DEG_progression_full.csv", "AAH vs LUAD (progression window)", "Fig2B_volcano_progression.png")

# ============ 2. 桥接基因热图 ============
d <- readRDS("step1_data_clean.rds")
expr <- d$expr; group <- d$group
bg <- bridge[bridge %in% rownames(expr)]
m <- expr[bg, , drop = FALSE]
z <- t(scale(t(m)))
z[z > 2] <- 2; z[z < -2] <- -2

ord <- order(group, colnames(z))
df <- data.frame(gene = factor(rep(rownames(z), ncol(z)), levels = rev(bg)),
                 sample = rep(colnames(z)[ord], each = nrow(z)),
                 value = as.vector(z[, ord]),
                 stage = factor(rep(as.character(group[ord]), each = nrow(z)),
                                levels = c("NL", "AAH", "LUAD")))
df$sample <- factor(df$sample, levels = colnames(z)[ord])

p2 <- ggplot(df, aes(sample, gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2E6DA4", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-2, 2), name = "z-score") +
  facet_grid(~ stage, scales = "free_x", space = "free_x") +
  labs(title = "Bridge genes (Rg3 targets that change across stages)",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), strip.background = element_rect(fill = "grey93"))
ggsave("Fig2C_bridge_heatmap.png", p2, width = 7.4, height = 4.6, dpi = 300)
cat("已保存 Fig2C_bridge_heatmap.png\n")

# ============ 3. PPI 网络（缺边文件就自动从 STRING 下载）============
if (!file.exists("PPI_string_edges.tsv")) {
  cat("未找到 PPI_string_edges.tsv，正在从 STRING 下载...\n")
  ids <- URLencode(paste(tg, collapse = "\r"), reserved = TRUE)
  url <- paste0("https://string-db.org/api/tsv/network?identifiers=", ids,
                "&species=9606&required_score=400&caller_identity=iris_rg3")
  ok <- tryCatch({ download.file(url, "PPI_string_edges.tsv", quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok || file.size("PPI_string_edges.tsv") < 100)
    stop("STRING 下载失败。请手动在浏览器打开下面这个链接，把页面内容另存为 PPI_string_edges.tsv：\n", url)
  cat("已下载：", file.size("PPI_string_edges.tsv"), "bytes\n")
}

e <- read.delim("PPI_string_edges.tsv")
e <- e[!duplicated(t(apply(e[, c("preferredName_A", "preferredName_B")], 1, sort))), ]
g <- graph_from_data_frame(e[, c("preferredName_A", "preferredName_B")], directed = FALSE)
deg <- sort(degree(g), decreasing = TRUE)
cat("PPI：", vcount(g), "个节点，", ecount(g), "条边\n")
is_bridge <- names(deg) %in% bg

V(g)$color <- ifelse(is_bridge, "#C0392B", "#9FB6C9")
V(g)$size  <- 2.5 + 0.40 * deg[names(V(g))]
V(g)$frame.color <- "white"
E(g)$width <- 0.3 + 0.8 * (e$score / max(e$score))

labA <- ifelse(is_bridge, names(deg), NA)
png("Fig2D_PPI_network.png", width = 3000, height = 1500, res = 260)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5))
set.seed(42)
plot(g, layout = layout_with_graphopt(g, niter = 800, charge = 0.02),
     vertex.label = labA, vertex.label.cex = 0.62, vertex.label.color = "black",
     vertex.label.family = "sans", edge.color = adjustcolor("#7F8C9B", 0.35),
     main = paste0("A  Full Rg3 target network (n = ", vcount(g), ")"))

core <- induced_subgraph(g, which(degree(g) >= 10))
cat("核心子网络：", vcount(core), "个节点，", ecount(core), "条边\n")
V(core)$size <- 3 + 0.5 * degree(core)
V(core)$color <- ifelse(names(V(core)) %in% bg, "#C0392B", "#9FB6C9")
plot(core, layout = layout_with_kk(core),
     vertex.label.cex = 0.85, vertex.label.color = "black",
     edge.color = adjustcolor("#7F8C9B", 0.5),
     main = paste0("B  Core subnetwork (degree >= 10, n = ", vcount(core), ")"))
legend("bottomleft", legend = c("Bridge gene (Rg3 target + differential)", "Other Rg3 target"),
       pt.bg = c("#C0392B", "#9FB6C9"), pch = 21, pt.cex = 1.8,
       bty = "n", cex = 0.8, title = "Node size = degree")
dev.off()
cat("已保存 Fig2D_PPI_network.png\n")

cat("\n>>> 四张图都生成了。\n")
