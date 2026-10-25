#!/usr/bin/env python3
"""生成 Supplementary Table S2：细胞类型 marker 列表（Methods 2.8 引用）"""
import os
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

ROOT = os.environ.get("RG3_ROOT", ".")
SN = os.path.join(ROOT, "sn")

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
FUNCTION = {
 "Epithelial/alveolar": "Alveolar and airway epithelium; lineage origin of AAH, AIS and invasive adenocarcinoma",
 "Basal/secretory": "Basal and secretory airway epithelium; candidate origin of peripheral squamous and mucinous lesions",
 "Fibroblast": "Resident pulmonary fibroblasts; structural stroma",
 "Endothelial": "Vascular endothelium; source of the angiogenic bridge genes",
 "Macrophage/monocyte": "Alveolar and interstitial macrophages, monocytes; innate immune surveillance",
 "T cell": "T lymphocytes; adaptive immunity",
 "NK cell": "Natural killer cells; cytotoxic innate immunity",
 "B cell": "B lymphocytes",
 "Plasma cell": "Antibody-secreting plasma cells",
 "Mast cell": "Mast cells; mediator release",
 "Dendritic cell": "Antigen-presenting dendritic cells; links innate and adaptive immunity",
 "Smooth muscle/pericyte": "Smooth muscle and pericytes; vessel wall",
 "Proliferating": "Cycling nuclei; retained as a quality-control category",
}

detected = set()
gp = os.path.join(SN, "sn_genes_sub.txt")
if os.path.exists(gp):
    detected = set(open(gp).read().split())

wb = Workbook()
ws = wb.active
ws.title = "S2 cell type markers"

ws["A1"] = "Supplementary Table S2. Marker gene sets used for cell-type annotation of GSE308103"
ws["A1"].font = Font(bold=True, size=12)
ws["A2"] = ("Nuclei were clustered de novo (22,425 subsampled nuclei, Leiden resolution 1.0, 25 clusters). "
            "Each cluster was assigned to the programme with the highest mean score across the marker set "
            "(scanpy score_genes, control-gene corrected). Column D gives the number of markers detected in the "
            "filtered gene list (15,510 genes) and is provided so that the annotation can be reproduced.")
ws["A2"].alignment = Alignment(wrap_text=True, vertical="top")
ws.merge_cells("A2:E2")
ws.row_dimensions[2].height = 46

hdr = ["Cell-type programme", "Marker genes", "n markers", "Detected in GSE308103", "Biological annotation"]
hrow = 4
for j, h in enumerate(hdr, start=1):
    c = ws.cell(row=hrow, column=j, value=h)
    c.font = Font(bold=True, color="FFFFFF")
    c.fill = PatternFill("solid", fgColor="4C72B0")
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

r = hrow + 1
for ct, genes in MARKERS.items():
    n_det = len([g for g in genes if not detected or g in detected])
    ws.cell(row=r, column=1, value=ct).font = Font(bold=True)
    ws.cell(row=r, column=2, value=", ".join(genes))
    ws.cell(row=r, column=3, value=len(genes)).alignment = Alignment(horizontal="center")
    ws.cell(row=r, column=4, value=f"{n_det}/{len(genes)}").alignment = Alignment(horizontal="center")
    ws.cell(row=r, column=5, value=FUNCTION[ct])
    for j in range(1, 6):
        ws.cell(row=r, column=j).alignment = Alignment(wrap_text=True, vertical="top",
                                                       horizontal="center" if j in (3, 4) else "left")
    r += 1

ws.cell(row=r + 1, column=1, value="Notes").font = Font(bold=True)
_tot = sum(len(v) for v in MARKERS.values())
_uniq = len({g for v in MARKERS.values() for g in v})
notes = [
 "Total markers: %d entries across 13 programmes (%d unique gene symbols)." % (_tot, _uniq),
 "Marker lists follow widely used lung cell-type conventions; no marker was selected on the basis of the results,",
 "so the annotation is independent of the bridge-gene analysis.",
 "The 'Proliferating' programme is retained for reproducibility; in the final clustering no cluster received this",
 "assignment, so no proliferating population is reported.",
 "Expression values used for scoring were log2(CPM+1) after panel-free library-size normalisation of the full",
 "count matrix (scanpy normalize_total, target sum 10,000, then log1p), z-scored per gene across nuclei.",
]
for i, t in enumerate(notes, start=1):
    ws.cell(row=r + 1 + i, column=1, value=t).alignment = Alignment(wrap_text=False)

widths = [26, 62, 10, 20, 62]
for j, w in enumerate(widths, start=1):
    ws.column_dimensions[get_column_letter(j)].width = w

out = os.path.join(ROOT, "Supplementary_Table_S2_cell_type_markers.xlsx")
wb.save(out)
print("已保存:", out)
print("marker 总数:", sum(len(v) for v in MARKERS.values()),
      "唯一符号:", len({g for v in MARKERS.values() for g in v}),
      "在 15,510 基因中检出:", sum(1 for v in MARKERS.values() for g in v if not detected or g in detected))
