#!/usr/bin/env python3
"""KEGG over-representation with an explicit, dataset-restricted background.

Background = KEGG-annotated protein-coding human genes ∩ genes retained in GSE102511.
Query sets = Rg3 predicted targets; the 22 bridge genes.
Test = two-sided Fisher exact test; BH-FDR within each query set.
Pathways with < 5 background genes are excluded.
"""
import csv, gzip, os, urllib.request
from math import comb

ROOT = os.environ.get("RG3_ROOT", ".")
CACHE = os.path.join(ROOT, "kegg_cache")
os.makedirs(CACHE, exist_ok=True)


def kegg(path):
    f = os.path.join(CACHE, path.replace("/", "_") + ".txt")
    if not os.path.exists(f):
        with urllib.request.urlopen("https://rest.kegg.jp/" + path, timeout=300) as r:
            open(f, "wb").write(r.read())
    return [l.rstrip("\n").split("\t") for l in open(f) if l.strip()]


names = {p: n for p, n in kegg("list/pathway/hsa")}
sym = {}  # hsa:ID -> official symbol (CDS only)
for row in kegg("list/hsa"):
    if len(row) >= 4 and row[1] == "CDS":
        s = row[3].split(";")[0].split(",")[0].strip()
        if s and " " not in s:
            sym[row[0]] = s
members = {}
for p, g in kegg("link/hsa/pathway"):
    p = p.replace("path:", "")
    if g in sym:
        members.setdefault(p, set()).add(sym[g])

# dataset genes = all genes tested by limma (DEG_onset_full.csv)
ds = {r["gene"] for r in csv.DictReader(open(os.path.join(ROOT, "DEG_onset_full.csv")))}
bg = set().union(*members.values()) & ds
print("dataset genes:", len(ds), "| KEGG background:", len(bg))
paths = {p: m & bg for p, m in members.items() if len(m & bg) >= 5}
print("testable pathways:", len(paths))


def fisher_two_sided(k, N, K, n):
    """Two-sided Fisher exact test on the 2x2 table (hit / pathway), as in scipy.stats.fisher_exact."""
    tot = comb(N, n)
    lo, hi = max(0, n - (N - K)), min(K, n)
    pmf = {i: comb(K, i) * comb(N - K, n - i) / tot for i in range(lo, hi + 1)}
    obs = pmf[k] * (1 + 1e-7)
    return min(1.0, sum(v for v in pmf.values() if v <= obs))


def bh(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q


def run(query, out):
    q = set(query) & bg
    N, n = len(bg), len(q)
    rows = []
    for p, m in paths.items():
        k = len(m & q)
        if k == 0:
            continue
        rows.append([p.replace("hsa", ""), names.get(p, p), k, len(m), fisher_two_sided(k, N, len(m), n)])
    fdr = bh([r[4] for r in rows])
    rows = [r + [f] for r, f in zip(rows, fdr)]
    rows.sort(key=lambda r: r[4])
    with open(os.path.join(ROOT, out), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pathway_id", "pathway", "n_hit", "n_pathway", "p", "fdr"]); w.writerows(rows)
    print(f"{out}: query {len(query)} -> mapped {n}; top = {rows[0][1]} (FDR {rows[0][5]:.2g})")


tg = {r["gene"] for r in csv.DictReader(open(os.path.join(ROOT, "rg3_targets_swisstarget.csv"))) if r["gene"] not in ("", "N/A", "N")}
br = [r["gene"] for r in csv.DictReader(open(os.path.join(ROOT, "bridge_genes.csv")))]
run(tg, "KEGG_targets93_bgDEG.csv")
run(br, "KEGG_bridge22_bgDEG.csv")
