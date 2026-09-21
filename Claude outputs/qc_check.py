"""
Verifica QC per Sergio - Visium 10x, campioni Skin/Wound1/Wound7/Wound30
(stesso donore, dataset GSE241124)

Due modalita' di calcolo QC richieste:
  A) per-sample: sc.pp.calculate_qc_metrics girato separatamente su ognuno
     dei 4 AnnData individuali -> plot separati per campione
  B) comprensivo: sc.pp.calculate_qc_metrics girato sull'AnnData unito
     (tutti gli spot insieme, 5955 x 36600) -> plot con confronto tra
     campioni (violin raggruppati) e plot con la distribuzione pooled
     (tutti gli spot mescolati, senza distinzione di campione)

Le metriche di gene-level (n_cells_by_counts, mean_counts) dipendono da QUALI
celle sono incluse nel calcolo: per questo motivo differiscono se calcolate
per-sample o sull'insieme -> utili per decidere una soglia di filtro sui geni
per-sample vs collettiva.
"""
import scanpy as sc
import spatialdata_io
import spatialdata as sd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

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


def add_mt_flag(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")


def qc_metrics(adata):
    add_mt_flag(adata)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata


# ---------------------------------------------------------------------
# A) caricamento + QC per-sample (calcolo indipendente su ogni campione)
# ---------------------------------------------------------------------
sdatas = {s: spatialdata_io.visium(f"{DATA}/{s}", dataset_id=s) for s in SAMPLES}
adatas = {}
for s in SAMPLES:
    a = sdatas[s]["table"].copy()
    qc_metrics(a)
    adatas[s] = a
    print(f"{s}: {a.n_obs} spot, {a.n_vars} geni | "
          f"median total_counts={a.obs['total_counts'].median():.0f}, "
          f"median n_genes={a.obs['n_genes_by_counts'].median():.0f}, "
          f"median pct_mt={a.obs['pct_counts_mt'].median():.2f}%")

# plot separato per campione (3 pannelli: total_counts, n_genes, pct_mt)
for s in SAMPLES:
    a = adatas[s]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, m in zip(axes, METRICS):
        sns.violinplot(y=a.obs[m], ax=ax, inner="box", cut=0)
        sns.stripplot(y=a.obs[m], ax=ax, color="black", size=1.5, alpha=0.3)
        ax.set_title(TITLES[m])
        ax.set_ylabel("")
    fig.suptitle(f"QC per-sample: {s}  (n={a.n_obs} spot)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/qc_violin_{s}.png", dpi=150)
    plt.close(fig)

# ---------------------------------------------------------------------
# B) QC comprensivo: calcolo unico sull'AnnData con tutti i campioni uniti
# ---------------------------------------------------------------------
sdata_all = sd.concatenate(list(sdatas.values()), concatenate_tables=True)
table_all = sdata_all["table"].copy()
qc_metrics(table_all)
table_all.obs["region"] = table_all.obs["region"].astype(str)

print()
print(f"Comprensivo (tutti insieme): {table_all.n_obs} spot totali, {table_all.n_vars} geni")
print(f"median total_counts={table_all.obs['total_counts'].median():.0f}, "
      f"median n_genes={table_all.obs['n_genes_by_counts'].median():.0f}, "
      f"median pct_mt={table_all.obs['pct_counts_mt'].median():.2f}%")

# B1) violin raggruppati per campione, calcolati UNA volta sull'oggetto unito
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for ax, m in zip(axes, METRICS):
    sns.violinplot(
        data=table_all.obs, x="region", y=m, order=SAMPLES,
        ax=ax, inner="box", cut=0,
    )
    ax.set_title(TITLES[m])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=30)
fig.suptitle("QC comprensivo: confronto tra campioni (calcolo unico su tutti gli spot)")
fig.tight_layout()
fig.savefig(f"{OUT}/qc_violin_combined_by_sample.png", dpi=150)
plt.close(fig)

# B2) distribuzione pooled: tutti gli spot mescolati, nessuna distinzione di campione
# (questa e' la vista da usare per decidere UNA soglia collettiva)
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, m in zip(axes, METRICS):
    sns.histplot(table_all.obs[m], bins=60, ax=ax, kde=True)
    ax.set_title(TITLES[m])
    ax.set_xlabel(m)
fig.suptitle(f"QC comprensivo: distribuzione pooled (n={table_all.n_obs} spot, tutti i campioni insieme)")
fig.tight_layout()
fig.savefig(f"{OUT}/qc_hist_combined_pooled.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# Gene-level: n_cells_by_counts / mean_counts, per-sample vs collettivo
# (serve a decidere la soglia min_cells per il filtro sui geni)
# ---------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for s in SAMPLES:
    sns.kdeplot(adatas[s].var["n_cells_by_counts"], ax=axes[0], label=s, log_scale=(True, False))
sns.kdeplot(table_all.var["n_cells_by_counts"], ax=axes[0], label="TUTTI (collettivo)",
            color="black", linewidth=2, linestyle="--", log_scale=(True, False))
axes[0].set_title("n_cells_by_counts per gene\n(in quanti spot e' rilevato ogni gene)")
axes[0].set_xlabel("n_cells_by_counts (log)")
axes[0].legend(fontsize=8)

for s in SAMPLES:
    sns.kdeplot(adatas[s].var["mean_counts"], ax=axes[1], label=s, log_scale=(True, False))
sns.kdeplot(table_all.var["mean_counts"], ax=axes[1], label="TUTTI (collettivo)",
            color="black", linewidth=2, linestyle="--", log_scale=(True, False))
axes[1].set_title("mean_counts per gene")
axes[1].set_xlabel("mean_counts (log)")
axes[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/qc_gene_metrics_per_sample_vs_combined.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# Salvataggio tabelle QC (obs/var) per ispezione numerica precisa
# ---------------------------------------------------------------------
obs_all = table_all.obs[["region", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
obs_all.to_csv(f"{OUT}/qc_obs_metrics_combined.csv")

var_all = table_all.var[["n_cells_by_counts", "mean_counts", "pct_dropout_by_counts"]].copy()
var_all.to_csv(f"{OUT}/qc_var_metrics_combined.csv")

summary_rows = []
for s in SAMPLES:
    o = adatas[s].obs
    summary_rows.append({
        "sample": s,
        "n_spots": adatas[s].n_obs,
        "median_total_counts": o["total_counts"].median(),
        "median_n_genes": o["n_genes_by_counts"].median(),
        "median_pct_mt": o["pct_counts_mt"].median(),
        "q10_total_counts": o["total_counts"].quantile(0.10),
        "q90_total_counts": o["total_counts"].quantile(0.90),
    })
summary = pd.DataFrame(summary_rows)
summary.to_csv(f"{OUT}/qc_summary_per_sample.csv", index=False)
print()
print(summary.to_string(index=False))

print()
print("File salvati in:", OUT)
import os
for f in sorted(os.listdir(OUT)):
    print(" -", f)
