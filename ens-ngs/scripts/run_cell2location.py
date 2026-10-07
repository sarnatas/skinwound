"""Cell2location deconvolution: reference signatures + spatial mapping.

Usage:
    python run_cell2location.py --profile vm  --threads 32 --inputs c2l_inputs --out c2l_out
    ppython run_cell2location.py --profile mac --threads 10 --inputs c2l_inputs --out c2l_out
"""
import argparse, os, time

parser = argparse.ArgumentParser()
parser.add_argument("--ref_epochs", type=int, default=None)
parser.add_argument("--vis_epochs", type=int, default=None)
parser.add_argument("--profile", choices=["vm", "mac"], required=True)
parser.add_argument("--threads", type=int, required=True)
parser.add_argument("--inputs", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

# Thread env vars MUST be set before importing torch / numpy-heavy libraries
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = str(args.threads)

import numpy as np
import pandas as pd
import anndata as ad
import torch
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel, Cell2location

torch.set_num_threads(args.threads)
scvi.settings.seed = 0
os.makedirs(args.out, exist_ok=True)

# ---- Profiles ---------------------------------------------------------------
PROFILES = {
    # full settings of the cell2location tutorial / the paper
    "vm":  dict(max_cells_per_state=None, ref_epochs=250, vis_epochs=30000, num_samples=1000),
    # quick version for a laptop: fewer cells, fewer epochs (check the ELBO curve!)
    "mac": dict(max_cells_per_state=1500, ref_epochs=150, vis_epochs=5000,  num_samples=500),
}
cfg = dict(PROFILES[args.profile])
if args.ref_epochs: cfg["ref_epochs"] = args.ref_epochs
if args.vis_epochs: cfg["vis_epochs"] = args.vis_epochs
ACC = "gpu" if torch.cuda.is_available() else "cpu"
print(f"profile={args.profile} accelerator={ACC} threads={args.threads} cfg={cfg}", flush=True)

# ---- 1. Reference -> cell-state signatures ---------------------------------------
t0 = time.time()
adata_ref = ad.read_h5ad(f"{args.inputs}/reference_counts.h5ad")

if cfg["max_cells_per_state"] is not None:
    keep = (adata_ref.obs.groupby("cell_state", observed=True, group_keys=False)
            .apply(lambda d: d.sample(min(len(d), cfg["max_cells_per_state"]), random_state=0))
            .index)
    adata_ref = adata_ref[adata_ref.obs_names.isin(keep)].copy()
print("reference:", adata_ref.shape, flush=True)

selected = filter_genes(adata_ref, cell_count_cutoff=5,
                        cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

RegressionModel.setup_anndata(adata=adata_ref, batch_key="sample", labels_key="cell_state")
mod_ref = RegressionModel(adata_ref)
mod_ref.train(max_epochs=cfg["ref_epochs"], accelerator=ACC)

fig = plt.figure(); mod_ref.plot_history(20); plt.savefig(f"{args.out}/ref_elbo.png", dpi=150); plt.close()

adata_ref = mod_ref.export_posterior(
    adata_ref, sample_kwargs={"num_samples": cfg["num_samples"], "batch_size": 2500, "accelerator": ACC})

factor_names = adata_ref.uns["mod"]["factor_names"]
inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{i}" for i in factor_names]].copy()
inf_aver.columns = factor_names
inf_aver.to_csv(f"{args.out}/reference_signatures.csv")
print(f"reference done in {(time.time()-t0)/60:.1f} min", flush=True)

# ---- 2. Spatial mapping -----------------------------------------------------------
t1 = time.time()
adata_vis = ad.read_h5ad(f"{args.inputs}/spatial_counts.h5ad")
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
print("spatial:", adata_vis.shape, flush=True)

Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")
mod_vis = Cell2location(adata_vis, cell_state_df=inf_aver,
                        N_cells_per_location=20, detection_alpha=20)
mod_vis.train(max_epochs=cfg["vis_epochs"], batch_size=None, train_size=1, accelerator=ACC)

fig = plt.figure(); mod_vis.plot_history(1000); plt.savefig(f"{args.out}/vis_elbo.png", dpi=150); plt.close()

adata_vis = mod_vis.export_posterior(
    adata_vis, sample_kwargs={"num_samples": cfg["num_samples"],
                              "batch_size": mod_vis.adata.n_obs, "accelerator": ACC})

# ---- 3. Save results ---------------------------------------------------------------
q05 = adata_vis.obsm["q05_cell_abundance_w_sf"].copy()
q05.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in q05.columns]
q05.to_csv(f"{args.out}/q05_cell_abundance.csv")           # small, easy to carry back
mod_vis.save(f"{args.out}/model_vis", overwrite=True)
adata_vis.write_h5ad(f"{args.out}/adata_vis_c2l.h5ad")
print(f"spatial done in {(time.time()-t1)/60:.1f} min", flush=True)