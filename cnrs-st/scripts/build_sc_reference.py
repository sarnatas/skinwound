"""Builds data/DATA_CNRS/sc_reference.h5ad (same logic as the notebook appendix) and checks the Visium barcodes."""
import sys, pandas as pd, anndata as ad, scanpy as sc
D = sys.argv[1]  # path to DATA_CNRS
if len(sys.argv) > 2 and sys.argv[2] == "--skip-check": meta = pd.DataFrame(columns=["Sample_name"])
# check: every barcode of the authors' Visium metadata is in the matching folder
meta = pd.read_csv(f"{D}/visium/GSE241124_spatialseq_metadata_acutewound.txt", sep="\t", index_col=0)
found = 0
for s, g in meta.groupby("Sample_name"):
    a = sc.read_10x_h5(f"{D}/visium/{s}/filtered_feature_bc_matrix.h5")
    ids = set(s + "_" + a.obs_names); found += g.index.isin(ids).sum()
print(f"Visium barcodes found: {found}/{len(meta)}", flush=True)
# scRNA reference (memory-friendly: choose the cells first, then read each sample and subset immediately)
sc_meta = pd.read_csv(f"{D}/scrna/GSE241132_cell_metadata.txt", sep="\t", index_col=0)
sc_meta = sc_meta[sc_meta["Patient"].isin(["PWH27", "PWH28"])]
keep = set(sc_meta.groupby("newCellTypes", group_keys=False).apply(lambda d: d.sample(min(len(d), 300), random_state=0)).index)
parts = []
for p in ["PWH27", "PWH28"]:
    for t in ["D0", "D1", "D7", "D30"]:
        a = sc.read_10x_mtx(f"{D}/scrna/{p}{t}"); a.var_names_make_unique()
        a.obs_names = [f"{p}{t}_{b}" for b in a.obs_names]
        parts.append(a[a.obs_names.isin(keep)].copy()); del a
ref = ad.concat(parts)
ref.obs = ref.obs.join(sc_meta[["Patient", "Condition", "newCellTypes", "newMainCellTypes"]])
ref.write_h5ad(f"{D}/sc_reference.h5ad", compression="gzip")
print(ref)
