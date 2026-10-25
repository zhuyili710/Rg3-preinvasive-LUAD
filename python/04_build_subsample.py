#!/usr/bin/env python3
"""第 1 步：把 75 个样本各抽 300 个核 -> 稀疏计数矩阵（int32，省内存）"""
import gzip, os, re, subprocess, tempfile, time
import numpy as np
import pandas as pd
import scipy.sparse as sp

SN = os.path.join(os.environ.get("RG3_ROOT", "."), "sn")
EX = os.path.join(SN, "extracted")
NCELL = 300

def parse_name(fn):
    m = re.match(r"GSM\d+_(P\d+)_([A-Za-z]+)", fn)
    return (m.group(1), re.sub(r"\d+$", "", m.group(2).upper())) if m else ("P?", "UNK")

files = sorted(f for f in os.listdir(EX) if f.endswith(".gz"))
print("样本数:", len(files), "每样本抽", NCELL, flush=True)

gene2i, chunks, cells, stages = {}, [], [], []
t0 = time.time()
for k, fn in enumerate(files):
    path = os.path.join(EX, fn)
    with gzip.open(path, "rt", errors="ignore") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
    stride = max(1, len(hdr) // NCELL)
    idx = list(range(0, len(hdr), stride))[:NCELL]
    fields = ",".join(["1"] + [str(i + 1) for i in idx])
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False); tmp.close()
    with open(tmp.name, "w") as out:
        subprocess.run(f'zcat "{path}" | cut -f{fields}', shell=True, stdout=out, check=True)
    df = pd.read_csv(tmp.name, sep="\t", header=None, index_col=0, skiprows=1, engine="c")
    os.unlink(tmp.name)
    df.index = df.index.astype(str)
    df = df[~df.index.duplicated()]
    for g in df.index:
        if g not in gene2i:
            gene2i[g] = len(gene2i)
    vals = np.nan_to_num(df.to_numpy(), nan=0.0).astype(np.int32)
    coo = sp.coo_matrix(vals)
    keep = coo.data > 0
    rows = coo.row[keep]; cols = coo.col[keep]; data = coo.data[keep]
    gidx = np.array([gene2i[df.index[r]] for r in rows], dtype=np.int32)
    chunks.append(sp.csr_matrix((data, (cols, gidx)), shape=(vals.shape[1], len(gene2i))))
    pid, stage = parse_name(fn)
    cells.extend([f"{pid}_{stage}_{k}_{j}" for j in range(vals.shape[1])])
    stages.extend([stage] * vals.shape[1])
    del df, vals, coo, rows, cols, data, gidx
    if k % 10 == 0 or k == len(files) - 1:
        print(f"  [{k+1}/{len(files)}] {pid} {stage} 基因={len(gene2i)} 用时={time.time()-t0:.0f}s", flush=True)

print("拼装中……", flush=True)
X = sp.vstack(chunks, format="csr")
del chunks
print("矩阵:", X.shape, "非零:", X.nnz, flush=True)

det = np.asarray((X > 0).sum(axis=0)).ravel()
keep_set = set(np.where(det >= 20)[0].tolist())
keep_g = sorted(keep_set)
X = X[:, keep_g].tocsr()
inv = {i: g for g, i in gene2i.items()}
genes = [inv[i] for i in keep_g]
print("过滤后:", X.shape, "非零:", X.nnz, flush=True)

sp.save_npz(os.path.join(SN, "sn_counts_sub.npz"), X)
with open(os.path.join(SN, "sn_genes_sub.txt"), "w") as f:
    f.write("\n".join(genes))
pd.DataFrame({"cell": cells, "stage": stages}).to_csv(os.path.join(SN, "sn_cells_sub.csv"), index=False)
print("完成，用时", round(time.time() - t0), "秒")
