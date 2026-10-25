#!/usr/bin/env python3
"""GSE308103 -> 每样本 pseudobulk（只累加目标基因，速度提升上百倍）

为什么要筛基因：这些 .mtx.txt 是"基因在行、细胞在列"的稠密文本，
全基因求和要解析上亿个数字，跑不完。我们只需要：
  - 22 个桥接基因
  - 两个窗口的差异基因
  - ME18 模块基因
其余基因对该节结论没有贡献。
"""
import gzip, os, re, csv
from collections import defaultdict

SN = os.path.join(os.environ.get("RG3_ROOT", "."), "sn")
EX = os.path.join(SN, "extracted")
ROOT = os.environ.get("RG3_ROOT", ".")

# ---- 目标基因表 ----
targets = set()
for f, col in [("bridge_genes.csv", "gene"), ("DEG_onset_sig.csv", "gene"),
               ("DEG_progression_sig.csv", "gene"), ("WGCNA_top_module_genes.csv", "gene")]:
    with open(os.path.join(ROOT, f), newline="") as fh:
        for r in csv.DictReader(fh):
            targets.add(r[col])
print("目标基因数:", len(targets), flush=True)

def parse_name(fn):
    base = os.path.basename(fn).replace(".raw_counts.mtx.txt.gz", "")
    m = re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)", base)
    if not m: return base, "UNK"
    return m.group(1), re.sub(r"\d+$", "", m.group(2).upper())

files = sorted(f for f in os.listdir(EX) if f.endswith(".gz"))
meta, pb = [], {}
for i, fn in enumerate(files):
    path = os.path.join(EX, fn)
    pid, stage = parse_name(fn)
    sums = defaultdict(int); ncells = 0
    with gzip.open(path, "rt", errors="ignore") as f:
        header = f.readline()
        ncells = header.count("\t")
        for line in f:
            t = line.find("\t")
            if t < 0: continue
            g = line[:t]
            if g not in targets: continue
            sums[g] += sum(map(int, line[t+1:].rstrip("\n").split("\t")))
    sample = fn.replace(".raw_counts.mtx.txt.gz", "")
    pb[sample] = sums
    meta.append({"sample": sample, "patient": pid,
                 "stage": stage, "cells": ncells})
    if i % 5 == 0 or i == len(files) - 1:
        print(f"  [{i+1}/{len(files)}] {pid} {stage} cells={ncells} hit={len(sums)}", flush=True)

genes = sorted(targets)
with open(os.path.join(SN, "sn_pseudobulk_selected.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene"] + [m["sample"] for m in meta])
    for g in genes:
        w.writerow([g] + [pb[m["sample"]].get(g, 0) for m in meta])
with open(os.path.join(SN, "sn_sample_meta.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["sample", "patient", "stage", "cells"])
    w.writeheader(); w.writerows(meta)
from collections import Counter
print("阶段分布:", Counter(m["stage"] for m in meta))
print("完成")
