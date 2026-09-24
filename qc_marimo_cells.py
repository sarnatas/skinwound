import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import scanpy as sc
    import spatialdata as sd
    import spatialdata_io
    from matplotlib import pyplot as plt
    import seaborn as sns
    import pandas as pd

    sns.set_style("whitegrid")
    return pd, plt, sc, sd, sns, spatialdata_io


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## QC — metriche per-sample e comprensive

    `sc.pp.calculate_qc_metrics` viene girato due volte: una volta per campione
    (per confrontare le distribuzioni tra Skin/Wound1/Wound7/Wound30), e una volta
    sull'oggetto con tutti gli spot uniti (per una vista comprensiva/pooled, utile
    per decidere soglie collettive).
    """)
    return


@app.cell
def _():
    data_folder = "/Users/sergio/Spatiotemporal_roadmap_of_human_skin_wound_healing/data/GSE241124_RAW"
    out_folder = "/Users/sergio/Spatiotemporal_roadmap_of_human_skin_wound_healing/qc_outputs"
    samples = ["Skin", "Wound1", "Wound7", "Wound30"]
    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    titles = {
        "total_counts": "UMI totali per spot",
        "n_genes_by_counts": "Geni rilevati per spot",
        "pct_counts_mt": "% conteggi mitocondriali per spot",
    }
    return data_folder, metrics, out_folder, samples, titles


@app.cell
def _(sc):
    def qc_metrics(adata):
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        return adata

    return (qc_metrics,)


@app.cell
def _(data_folder, qc_metrics, samples, spatialdata_io):
    # QC per-sample: calcolo indipendente su ciascun campione
    adatas = {}
    for _s in samples:
        _sdata = spatialdata_io.visium(f"{data_folder}/{_s}", dataset_id=_s)
        _a = _sdata["table"].copy()
        qc_metrics(_a)
        adatas[_s] = _a
    return (adatas,)


@app.cell
def _(data_folder, qc_metrics, samples, sd, spatialdata_io):
    # QC comprensivo: calcolo unico sull'oggetto con tutti gli spot uniti
    _sdatas = {
        _s: spatialdata_io.visium(f"{data_folder}/{_s}", dataset_id=_s) for _s in samples
    }
    sdata_all = sd.concatenate(list(_sdatas.values()), concatenate_tables=True)
    table_all = sdata_all["table"].copy()
    qc_metrics(table_all)
    table_all.obs["region"] = table_all.obs["region"].astype(str)
    return (table_all,)


@app.cell
def _(adatas, metrics, out_folder, plt, samples, sns, titles):
    # plot A: violin per-sample, griglia campioni (righe) x metriche (colonne)
    # tutte le variabili di lavoro con underscore -> riusabili nelle celle sotto
    _fig, _axes = plt.subplots(
        len(samples), len(metrics), figsize=(12, 3.2 * len(samples))
    )
    for _row, _s in enumerate(samples):
        _a = adatas[_s]
        for _col, _m in enumerate(metrics):
            _ax = _axes[_row, _col]
            sns.violinplot(y=_a.obs[_m], ax=_ax, inner="box", cut=0)
            sns.stripplot(y=_a.obs[_m], ax=_ax, color="black", size=1.5, alpha=0.3)
            _ax.set_ylabel("")
            if _row == 0:
                _ax.set_title(titles[_m])
            if _col == 0:
                _ax.set_ylabel(f"{_s}\n(n={_a.n_obs})", fontsize=11, fontweight="bold")
    _fig.suptitle("QC per-sample: confronto verticale (un calcolo per campione)", y=1.0)
    _fig.tight_layout()
    _fig.savefig(f"{out_folder}/qc_violin_all_samples_stacked.png", dpi=150, bbox_inches="tight")
    _fig
    return


@app.cell
def _(metrics, out_folder, plt, samples, sns, table_all, titles):
    # plot B: violin comprensivo, un solo calcolo, raggruppato per campione
    _fig, _axes = plt.subplots(1, len(metrics), figsize=(14, 4.5))
    for _ax, _m in zip(_axes, metrics):
        sns.violinplot(
            data=table_all.obs, x="region", y=_m, order=samples, ax=_ax, inner="box", cut=0
        )
        _ax.set_title(titles[_m])
        _ax.set_xlabel("")
        _ax.set_ylabel("")
        _ax.tick_params(axis="x", rotation=30)
    _fig.suptitle("QC comprensivo: confronto tra campioni (calcolo unico su tutti gli spot)")
    _fig.tight_layout()
    _fig.savefig(f"{out_folder}/qc_violin_combined_by_sample.png", dpi=150)
    _fig
    return


@app.cell
def _(metrics, out_folder, plt, sns, table_all):
    # plot C: distribuzione pooled, tutti gli spot mescolati, nessuna distinzione di campione
    # -> questa è la vista da guardare per decidere UNA soglia collettiva
    _fig, _axes = plt.subplots(1, len(metrics), figsize=(14, 4))
    for _ax, _m in zip(_axes, metrics):
        sns.histplot(table_all.obs[_m], bins=60, ax=_ax, kde=True)
        _ax.set_title(_m)
        _ax.set_xlabel(_m)
    _fig.suptitle(f"QC comprensivo: distribuzione pooled (n={table_all.n_obs} spot)")
    _fig.tight_layout()
    _fig.savefig(f"{out_folder}/qc_hist_combined_pooled.png", dpi=150)
    _fig
    return


@app.cell
def _(adatas, out_folder, plt, samples, sns, table_all):
    # plot D: gene-level, n_cells_by_counts / mean_counts, per-sample vs collettivo
    # -> utile per la soglia min_cells del filtro sui geni
    _fig, _axes = plt.subplots(1, 2, figsize=(11, 4.5))
    _palette = sns.color_palette("tab10", n_colors=len(samples))

    for _s, _c in zip(samples, _palette):
        _vals = adatas[_s].var["n_cells_by_counts"]
        _vals = _vals[_vals > 0]
        sns.histplot(
            _vals, ax=_axes[0], color=_c, label=_s, stat="density",
            element="step", fill=False, log_scale=True, bins=40,
        )
    _vals_all = table_all.var["n_cells_by_counts"]
    _vals_all = _vals_all[_vals_all > 0]
    sns.histplot(
        _vals_all, ax=_axes[0], color="black", label="TUTTI (collettivo)", stat="density",
        element="step", fill=False, log_scale=True, bins=40, linewidth=2, linestyle="--",
    )
    _axes[0].set_title("n_cells_by_counts per gene (>0)")
    _axes[0].legend(fontsize=8)

    for _s, _c in zip(samples, _palette):
        _vals = adatas[_s].var["mean_counts"]
        _vals = _vals[_vals > 0]
        sns.histplot(
            _vals, ax=_axes[1], color=_c, label=_s, stat="density",
            element="step", fill=False, log_scale=True, bins=40,
        )
    _vals_all = table_all.var["mean_counts"]
    _vals_all = _vals_all[_vals_all > 0]
    sns.histplot(
        _vals_all, ax=_axes[1], color="black", label="TUTTI (collettivo)", stat="density",
        element="step", fill=False, log_scale=True, bins=40, linewidth=2, linestyle="--",
    )
    _axes[1].set_title("mean_counts per gene (>0)")
    _axes[1].legend(fontsize=8)

    _fig.tight_layout()
    _fig.savefig(f"{out_folder}/qc_gene_metrics_per_sample_vs_combined.png", dpi=150)
    _fig
    return


@app.cell
def _(adatas, out_folder, pd, samples):
    # tabella riassuntiva per-sample, salvata come CSV
    _rows = []
    for _s in samples:
        _o = adatas[_s].obs
        _rows.append(
            {
                "sample": _s,
                "n_spots": adatas[_s].n_obs,
                "median_total_counts": _o["total_counts"].median(),
                "median_n_genes": _o["n_genes_by_counts"].median(),
                "median_pct_mt": _o["pct_counts_mt"].median(),
            }
        )
    summary = pd.DataFrame(_rows)
    summary.to_csv(f"{out_folder}/qc_summary_per_sample.csv", index=False)
    summary
    return


if __name__ == "__main__":
    app.run()
