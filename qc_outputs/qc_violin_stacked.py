import scanpy as sc
import spatialdata_io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DATA = "/mnt/user-data/uploads/Spatiotemporal_roadmap_of_human_skin_wound_healing/data/GSE241124_RAW"
OUT = "/tmp/claude-0/-home-claude/d839a57e-d5ba-5fda-95b3-059d2a571b65/scratchpad/qc_outputs"
SAMPLES = ["Skin", "Wound1", "Wound7", "Wound30"]
METRICS = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
TITLES = {
    "total_counts": "UMI totali per spot",
    "n_genes_by_counts": "Geni rilevati per spot",
    "pct_counts_mt": "% conteggi mitocondriali per spot",
}
sns.set_style("whitegrid")

def qc_metrics(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata

adatas = {}
for s in SAMPLES:
    sdata = spatialdata_io.visium(f"{DATA}/{s}", dataset_id=s)
    a = sdata["table"].copy()
    qc_metrics(a)
    adatas[s] = a

# --- griglia SAMPLES (righe) x METRICS (colonne), un pannello per cella ---
fig, axes = plt.subplots(len(SAMPLES), len(METRICS), figsize=(12, 3.2 * len(SAMPLES)))

for row, s in enumerate(SAMPLES):
    a = adatas[s]
    for col, m in enumerate(METRICS):
        ax = axes[row, col]
        sns.violinplot(y=a.obs[m], ax=ax, inner="box", cut=0)
        sns.stripplot(y=a.obs[m], ax=ax, color="black", size=1.5, alpha=0.3)
        ax.set_ylabel("")
        if row == 0:
            ax.set_title(TITLES[m])
        if col == 0:
            ax.set_ylabel(f"{s}\n(n={a.n_obs})", fontsize=11, fontweight="bold")

fig.suptitle("QC per-sample: Skin / Wound1 / Wound7 / Wound30 (un calcolo per campione)", y=1.0)
fig.tight_layout()
fig.savefig(f"{OUT}/qc_violin_all_samples_stacked.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Salvato:", f"{OUT}/qc_violin_all_samples_stacked.png")
