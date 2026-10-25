# Ginsenoside Rg3 targets the pre-invasive transition of lung adenocarcinoma

Companion code and intermediate result tables for the manuscript:

> Ginsenoside Rg3 targets the pre-invasive transition of lung adenocarcinoma: an integrative
> analysis of bulk and single-nucleus transcriptomes with network pharmacology.

## Data (all public; none redistributed here except small derived tables)

| Dataset | Use | Access |
|---|---|---|
| GSE102511 | Bulk RNA-seq, 17 donors x 3 stages (normal / AAH / LUAD) | downloaded automatically by `R/01_load_qc.R` |
| GSE308103 | snRNA-seq validation (normal / AAH / AIS / MIA / LUAD) | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE308nnn/GSE308103/suppl/GSE308103_RAW.tar` (1.5 GB) |
| SwissTargetPrediction | Rg3 (PubChem CID 9918693) predicted targets | `rg3_targets_swisstarget.csv` (included, repo root) |
| STRING v12 | protein-protein interaction edges, score >= 400 | `PPI_string_edges.tsv` (included, repo root) / re-downloaded by `R/04_figures.R` |
| KEGG REST API | pathway definitions | `rest.kegg.jp`, cached by `python/01_kegg_enrichment.py` |

## Order of execution

```
export RG3_ROOT=/path/to/this/repo

# 1. Bulk RNA-seq: QC, paired limma DE, WGCNA, figures
Rscript R/01_load_qc.R
Rscript R/02_limma_two_windows.R      # -> DEG_onset_*.csv, DEG_progression_*.csv
Rscript R/03_wgcna.R                  # -> WGCNA_*.csv
Rscript R/04_figures.R                # -> Fig2A-D; computes the 22 bridge genes

# 2. KEGG over-representation (two-sided Fisher exact, dataset-restricted background)
python3 python/01_kegg_enrichment.py  # -> KEGG_targets93_bgDEG.csv, KEGG_bridge22_bgDEG.csv

# 3. Single-nucleus validation (GSE308103)
#    download GSE308103_RAW.tar into $RG3_ROOT/sn/, then: tar xf GSE308103_RAW.tar -C sn/extracted
python3 python/02_pseudobulk.py        # per-sample pseudobulk, 999 target genes only (speed)
python3 python/03_validate_direction.py  # direction concordance; Spearman rho
Rscript R/05_fig3_snRNA.R              # -> Fig3_snRNA_validation.png

# 4. Cell-type localisation (22,425-nucleus subsample; full series exceeded available memory)
python3 python/04_build_subsample.py   # -> sn_counts_sub.npz (300 nuclei/sample, fixed stride)
python3 python/05_cluster_annotate.py  # Leiden clustering + marker-score cell-type calls

# 5. Supplementary Table S2
python3 python/06_make_table_S2.py
```

## Reference outputs

`results/` contains the exact tables used in the manuscript (DEG lists, bridge genes, WGCNA module–stage
correlations, PPI degree, KEGG enrichment, single-nucleus composition and bridge-gene expression by
cell type and stage). Re-running the pipeline should reproduce them; KEGG output was verified to be
identical when re-run from this repository.

## Environment

- R 4.3.2+ with limma 3.58.1, WGCNA 1.72.5, ggplot2, igraph, ggrepel, data.table
- Python 3.12 with scanpy 1.12.4, anndata, leidenalg, python-igraph, umap-learn, pandas, numpy, scipy, matplotlib, openpyxl
- STRING and KEGG were accessed 22 October 2026; pathway/interaction content may change on re-query

## Notes on analysis choices

- Differential expression uses a donor-paired design (`~0 + stage + donor`); group contrasts are
  built with `makeContrasts()`. Passing a coefficient name directly to `topTable()` on a cell-means
  design returns absolute expression, not a group difference — this is a common limma pitfall.
- WGCNA runs single-threaded (`disableWGCNAThreads()`); on a memory-constrained machine,
  `maxBlockSize` was capped at 5000 and soft-threshold selection used only the 5000 most variable
  genes to avoid out-of-memory kills.
- Gene symbols in the public GSE102511 matrix include spreadsheet-mangled entries
  (e.g. `MARCH1` stored as `"1-Mar"`); these are repaired in `R/03_wgcna.R`.
- KEGG enrichment uses an explicit background (KEGG-annotated genes intersected with the genes
  retained after QC, not all KEGG genes) and a two-sided Fisher exact test, matching the rest of
  the manuscript's statistical convention. (A one-sided hypergeometric
  test would give smaller p-values for large pathways; we report the two-sided form.)
- The single-nucleus cell-type analysis is a 22,425-nucleus subsample (300 nuclei per sample,
  fixed stride), not the full 798,025-nucleus series, for the same memory reason as WGCNA.
  Cluster proportions should be read as estimates, not exact counts.

## License

Code: MIT. Data tables derive entirely from public sources (GEO, STRING, KEGG, SwissTargetPrediction)
and are redistributed here only as the small derived summaries needed to reproduce the figures.
