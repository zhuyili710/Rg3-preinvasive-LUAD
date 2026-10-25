#!/usr/bin/env python3
"""3.8 单核验证：GSE308103 pseudobulk 层面的阶段一致性检验"""
import csv, os
from collections import defaultdict
import statistics as st

SN = os.path.join(os.environ.get("RG3_ROOT", "."), "sn")
ROOT = os.environ.get("RG3_ROOT", ".")
ORDER = ["NORMAL", "AAH", "AIS", "MIA", "LUAD"]

meta = list(csv.DictReader(open(os.path.join(SN, "sn_sample_meta.csv"))))
rows = list(csv.reader(open(os.path.join(SN, "sn_pseudobulk_selected.csv"))))
genes = [r[0] for r in rows[1:]]
mat = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
samples = rows[0][1:]
print("样本", len(samples), "基因", len(genes))
from collections import Counter
print("阶段分布:", Counter(m["stage"] for m in meta))

# 面板内 CPM（用 999 个分析基因的合计做文库大小代理）
tot = [sum(mat[g][i] for g in genes) for i in range(len(samples))]
logn = {g: [__import__("math").log2(mat[g][i] / tot[i] * 1e6 + 1) for i in range(len(samples))]
        for g in genes}

def mean_by_stage(v):
    d = defaultdict(list)
    for i, m in enumerate(meta):
        d[m["stage"]].append(v[i])
    return {k: st.mean(d[k]) for k in ORDER if k in d}

def spearman(x, y):
    def rank(a):
        idx = sorted(range(len(a)), key=lambda i: a[i]); r = [0]*len(a)
        i = 0
        while i < len(a):
            j = i
            while j+1 < len(a) and a[idx[j+1]] == a[idx[i]]: j += 1
            for k in range(i, j+1): r[idx[k]] = (i+j)/2 + 1
            i = j+1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x); mx, my = st.mean(rx), st.mean(ry)
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = (sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry)) ** 0.5
    return num/den if den else 0.0

# ---- 1) ME18（阶段负相关模块）在单核 pseudobulk 中是否同样下降 ----
me18 = [r["gene"] for r in csv.DictReader(open(os.path.join(ROOT, "WGCNA_top_module_genes.csv")))]
me18 = [g for g in me18 if g in logn]
sig = [st.mean([logn[g][i] for g in me18]) for i in range(len(samples))]
mbs = mean_by_stage(sig)
stage_rank = [ORDER.index(m["stage"]) for m in meta if m["stage"] in ORDER]
keep = [i for i, m in enumerate(meta) if m["stage"] in ORDER]
rho = spearman([stage_rank[k] for k in range(len(stage_rank))], [sig[i] for i in keep])
print(f"\n[ME18] {len(me18)} 个基因 → 各阶段均值:", {k: round(v, 2) for k, v in mbs.items()})
print(f"[ME18] 与阶段秩的 Spearman rho = {rho:.3f}")

# ---- 2) 全局方向一致性：bulk logFC vs 单核阶段差 ----
def bulk(f):
    return {r["gene"]: float(r["logFC"]) for r in csv.DictReader(open(os.path.join(ROOT, f)))}

for (g1, g2, bf, label) in [("NORMAL", "AAH", "DEG_onset_full.csv", "onset  NORMAL->AAH"),
                            ("AAH", "LUAD", "DEG_progression_full.csv", "prog   AAH->LUAD")]:
    B = bulk(bf)
    gs = [g for g in genes if g in B]
    xs, ys = [], []
    for g in gs:
        i1 = [i for i, m in enumerate(meta) if m["stage"] == g1]
        i2 = [i for i, m in enumerate(meta) if m["stage"] == g2]
        if not i1 or not i2: continue
        xs.append(B[g]); ys.append(st.mean([logn[g][i] for i in i2]) - st.mean([logn[g][i] for i in i1]))
    r = spearman(xs, ys)
    # 显著差异基因的方向一致率
    sig_genes = [g for g in gs if abs(B[g]) > 1]
    agree = 0
    for g in sig_genes:
        i1 = [i for i, m in enumerate(meta) if m["stage"] == g1]
        i2 = [i for i, m in enumerate(meta) if m["stage"] == g2]
        d = st.mean([logn[g][i] for i in i2]) - st.mean([logn[g][i] for i in i1])
        if B[g]*d > 0: agree += 1
    print(f"\n[{label}] 可核对基因 {len(xs)} 个；Spearman rho = {r:.3f}")
    print(f"[{label}] |logFC|>1 的 {len(sig_genes)} 个基因中方向一致：{agree} ({100*agree/max(1,len(sig_genes)):.1f}%)")

# ---- 3) 22 个桥接基因逐一核对 ----
bridge = [r["gene"] for r in csv.DictReader(open(os.path.join(ROOT, "bridge_genes.csv")))]
det = {r["gene"]: r for r in csv.DictReader(open(os.path.join(ROOT, "bridge_genes_detailed.csv")))}
print("\n[桥接基因] 单核 pseudobulk 阶段趋势：")
consist = 0
for g in bridge:
    if g not in logn: continue
    ms = mean_by_stage(logn[g])
    # 用该基因所属窗口的 bulk 方向做对照
    w = det.get(g, {}).get("window", "onset")
    b = bulk("DEG_onset_full.csv" if w == "onset" else "DEG_progression_full.csv")[g]
    d = (ms.get("AAH", ms.get("NORMAL")) - ms.get("NORMAL", 0)) if w == "onset" else (ms.get("LUAD", 0) - ms.get("AAH", 0))
    ok = (b*d) > 0
    consist += ok
    print(f"   {g:9s} {w:11s} bulk_logFC={b:+.2f} sn_delta={d:+.2f} {'一致' if ok else '不一致'}")
print(f"\n[桥接基因] 方向一致 {consist}/{len(bridge)}")
