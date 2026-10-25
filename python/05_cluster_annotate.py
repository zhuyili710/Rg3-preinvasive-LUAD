#!/usr/bin/env python3
"""第 2 步：Leiden 聚类 + marker 注释 + 桥接基因定位（抽样 2.25 万核）"""
import os, json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

sc.settings.verbosity = 1
SN = os.path.join(os.environ.get("RG3_ROOT", "."), "sn")
ROOT = os.environ.get("RG3_ROOT", ".")
OUT = os.path.join(SN, "cluster_out")
os.makedirs(OUT, exist_ok=True)

X = sp.load_npz(os.path.join(SN, "sn_counts_sub.npz"))
genes = open(os.path.join(SN, "sn_genes_sub.txt")).read().split()
cells = pd.read_csv(os.path.join(SN, "sn_cells_sub.csv"))
ad = sc.AnnData(X=X, obs=cells.set_index("cell"), var=pd.DataFrame(index=genes))
ad.obs["stage"] = pd.Categorical(ad.obs["stage"], categories=["NORMAL","AAH","AIS","MIA","LUAD"])
print("读入:", ad.shape, flush=True)

ad.layers["counts"] = ad.X.copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
print("HVG:", int(ad.var.highly_variable.sum()), flush=True)

hv = ad[:, ad.var.highly_variable].copy()
sc.pp.scale(hv, max_value=10)
sc.tl.pca(hv, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(hv, n_neighbors=15, n_pcs=20)
sc.tl.leiden(hv, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
ad.obs["leiden"] = hv.obs["leiden"].values
print("簇数:", ad.obs["leiden"].nunique(), flush=True)

try:
    sc.tl.umap(hv, min_dist=0.3)
    ad.obsm["X_umap"] = hv.obsm["X_umap"]
except Exception as e:
    print("UMAP 失败（不影响结论）:", e, flush=True)

MARKERS = {
 "Epithelial/alveolar": ["EPCAM","KRT8","KRT18","KRT19","SFTPC","SFTPB","SFTPA1","SFTPA2","AGER","NAPSA","CLDN18","MUC1"],
 "Basal/secretory":     ["KRT5","KRT14","TP63","KRT17","MUC5B","SCGB1A1","SCGB3A2","FOXJ1"],
 "Fibroblast":          ["COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","FBLN1"],
 "Endothelial":         ["PECAM1","VWF","CLDN5","EMCN","KDR","PLVAP"],
 "Macrophage/monocyte": ["CD68","CD14","MRC1","C1QA","C1QB","APOE","LYZ","ITGAM","MARCO","FABP4","PPARG"],
 "T cell":              ["CD3D","CD3E","CD2","CD8A","CD4","IL7R","TRAC","TRBC2"],
 "NK cell":             ["NKG7","GNLY","KLRD1","PRF1"],
 "B cell":              ["MS4A1","CD79A","CD79B","CD19","IGHM"],
 "Plasma cell":         ["MZB1","JCHAIN","IGHG1","SDC1"],
 "Mast cell":           ["TPSAB1","CPA3","KIT","MS4A2"],
 "Dendritic cell":      ["CD1C","CLEC9A","LAMP3","ITGAX"],
 "Smooth muscle/pericyte": ["ACTA2","RGS5","MYH11","TAGLN"],
 "Proliferating":       ["MKI67","TOP2A"],
}
import re as _re
safe = {ct: "sc_" + _re.sub(r"[^A-Za-z0-9_]", "_", ct) for ct in MARKERS}
for ct, gl in MARKERS.items():
    g = [x for x in gl if x in ad.var_names]
    if g:
        sc.tl.score_genes(ad, g, score_name=safe[ct], use_raw=False)

inv = {v: k for k, v in safe.items()}
scores = ad.obs[[c for c in ad.obs.columns if c.startswith("sc_")]]
cl_mean = scores.groupby(ad.obs["leiden"].values).mean()
assign = {cl: inv[cl_mean.columns[i]] for cl, i in
          zip(cl_mean.index, cl_mean.values.argmax(axis=1))}
ad.obs["celltype"] = ad.obs["leiden"].map(assign)

cl_mean.rename(columns=inv).round(2).to_csv(os.path.join(OUT, "cluster_marker_scores.csv"))
cols = [c for c in ["leiden","celltype","stage","sample"] if c in ad.obs.columns]
ad.obs[cols].to_csv(os.path.join(OUT, "cell_calls.csv"))

# 每个 cluster 的组成 + 每阶段比例
tab = pd.crosstab(ad.obs["celltype"], ad.obs["stage"], normalize="columns") * 100
tab.round(1).to_csv(os.path.join(OUT, "composition_by_stage.csv"))
print("\n各阶段细胞类型占比(%):\n", tab.round(1).to_string(), flush=True)

# 桥接基因按细胞类型 / 阶段
bridge = [g for g in pd.read_csv(os.path.join(ROOT, "bridge_genes.csv"))["gene"] if g in ad.var_names]
sub = ad[:, bridge]
df = pd.DataFrame(np.asarray(sub.X.todense()), columns=bridge)
df["celltype"] = ad.obs["celltype"].values
df["stage"] = ad.obs["stage"].values
df.groupby("celltype")[bridge].mean().round(3).to_csv(os.path.join(OUT, "bridge_by_celltype.csv"))
df.groupby(["celltype","stage"])[bridge].mean().round(3).to_csv(os.path.join(OUT, "bridge_by_celltype_stage.csv"))

print("\n桥接基因主要细胞类型:", flush=True)
m = df.groupby("celltype")[bridge].mean()
for g in bridge:
    print(f"  {g:9s} -> {m[g].idxmax():22s} {m[g].max():.2f}", flush=True)

# 图
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
if "X_umap" in ad.obsm:
    for key, fn in [("leiden","Fig4A_umap_clusters.png"), ("celltype","Fig4B_umap_celltype.png"), ("stage","Fig4C_umap_stage.png")]:
        sc.pl.umap(ad, color=key, legend_loc="on data" if key != "leiden" else "right margin",
                   show=False, frameon=False, title=key)
        plt.savefig(os.path.join(OUT, fn), dpi=250, bbox_inches="tight"); plt.close()
    sc.pl.dotplot(ad, bridge, groupby="celltype", standard_scale="var", show=False)
    plt.savefig(os.path.join(OUT, "Fig4D_bridge_dotplot.png"), dpi=250, bbox_inches="tight"); plt.close()
print("已保存图到", OUT, flush=True)
try:
    ad.obs["celltype"] = ad.obs["celltype"].astype(str)
    ad.obs["leiden"] = ad.obs["leiden"].astype(str)
    ad.write(os.path.join(OUT, "adata_clustered.h5ad"))
    print("已保存 adata_clustered.h5ad", flush=True)
except Exception as e:
    print("h5ad 保存失败（不影响其它结果）:", e, flush=True)
