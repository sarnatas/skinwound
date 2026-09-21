import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import numpy as np
    import pandas as pd

    return


@app.cell
def _():
    from matplotlib import pyplot as plt

    return (plt,)


@app.cell
def _():
    import scanpy as sc
    import squidpy as sq
    import spatialdata as sd
    import spatialdata_io
    import spatialdata_plot
    import seaborn as sns

    return sc, sns, spatialdata_io


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We will analyze data coming from one of the four donors reported in the paper.
    All the data will be loaded in the same spatialdata object and analyzed
    """)
    return


@app.cell
def _(sns):
    DATA = '/Users/sergio/Spatiotemporal_roadmap_of_human_skin_wound_healing/data/GSE241124_RAW'
    OUT = '/Users/sergio/Spatiotemporal_roadmap_of_human_skin_wound_healing/output/'
    SAMPLES = ["Skin", "Wound1", "Wound7", "Wound30"]
    METRICS = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    TITLES = {
        "total_counts": "Total UMIs per spot",
        "n_genes_by_counts": "Detected genes per spot",
        "pct_counts_mt": "% mitochondrial counts per spot",
    }

    sns.set_style("whitegrid")
    return DATA, METRICS, OUT, SAMPLES, TITLES


@app.cell
def _(
    DATA,
    METRICS,
    OUT,
    SAMPLES,
    TITLES,
    plt,
    qc_metrics,
    sns,
    spatialdata_io,
):
    # QC per-sample (independent for each sample)

    sdatas = {s: spatialdata_io.visium(f"{DATA}/{s}", dataset_id=s) for s in SAMPLES}
    adatas = {}

    for s in SAMPLES:
        a = sdatas[s]["table"].copy()
        qc_metrics(a)
        adatas[s] = a
        print(f"{s}: {a.n_obs} spot, {a.n_vars} geni | "
              f"mean total_counts={a.obs['total_counts'].mean():.0f}, "
              f"mean n_genes={a.obs['n_genes_by_counts'].mean():.0f}, "
              f"mean pct_mt={a.obs['pct_counts_mt'].mean():.2f}%")

    # Plot and save per sample
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
    return (sdatas,)


@app.cell
def _(samples, sdatas):
    # Check the common genes among the samples:
    gene_sets = {s: set(sdatas[s]["table"].var_names) for s in samples}
    common = set.intersection(*gene_sets.values())
    union = set.union(*gene_sets.values())
    print(len(common), len(union))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Quality control
    """)
    return


@app.cell
def _(sc):
    def add_mt_flag(adata):
        adata.var["mt"] = adata.var_names.str.startswith("MT-") 

    def qc_metrics(adata):
        add_mt_flag(adata)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        return adata


    return (qc_metrics,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
