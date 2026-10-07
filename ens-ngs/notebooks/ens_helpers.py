"""
ens_helpers.py - helper functions for the ENS-NGS practical
(spatial transcriptomics + scRNA-seq of human skin wound healing).

These functions wrap *plumbing* code (plots repeated for each sample, SpatialData
bookkeeping, small pandas reshaping) so that the student notebook can focus on the
biology and on the analysis choices. Nothing in here is "magic": every function has a
docstring that says exactly what it does. If you are curious, read the source!

Usage in the notebook:
    import ens_helpers as eh
    help(eh.qc_metrics)
"""

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from matplotlib.colors import Normalize


# ----------------------------------------------------------------------------
# Quality control
# ----------------------------------------------------------------------------
def qc_metrics(adata):
    """Compute the standard QC metrics, in place, and return the object.

    1. flags mitochondrial genes (gene symbol starting with "MT-") in adata.var["mt"]
    2. calls sc.pp.calculate_qc_metrics, which adds to adata.obs:
       total_counts, n_genes_by_counts, pct_counts_mt (and total_counts_mt)
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def plot_qc_violins(adatas, metrics, titles, unit="spots"):
    """One figure per sample: a violin + strip plot for each QC metric.

    adatas  : dict {sample name: AnnData}
    metrics : list of columns of .obs to plot
    titles  : dict {metric: panel title}
    """
    for name, a in adatas.items():
        fig, axes = plt.subplots(1, len(metrics), figsize=(12, 4))
        for ax, metric in zip(axes, metrics):
            sns.violinplot(y=a.obs[metric], ax=ax, inner="box", cut=0)
            sns.stripplot(y=a.obs[metric], ax=ax, color="black", size=1.5, alpha=0.3)
            ax.set_title(titles[metric])
            ax.set_ylabel("")
        fig.suptitle(f"QC per sample: {name} (n={a.n_obs} {unit})")
        fig.tight_layout()
        plt.show()


def plot_qc_combined(obs, group_col, order, metrics, titles, suptitle="QC combined across samples"):
    """One figure, one panel per QC metric, one violin per sample (side by side).

    obs       : a DataFrame (typically adata.obs)
    group_col : column of `obs` that contains the sample name
    order     : order of the samples on the x axis
    """
    fig, axes = plt.subplots(1, len(metrics), figsize=(14, 4))
    for ax, metric in zip(axes, metrics):
        sns.violinplot(data=obs, x=group_col, y=metric, order=order, ax=ax, inner="box", cut=0)
        ax.set_title(titles[metric])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=30)
    fig.suptitle(suptitle)
    fig.tight_layout()
    plt.show()


# ----------------------------------------------------------------------------
# SpatialData / squidpy bookkeeping
# ----------------------------------------------------------------------------
def attach_spatial_uns(adata, sdata, sample, data_folder):
    """Rebuild the classic adata.uns["spatial"] structure that squidpy expects.

    spatialdata_io.visium() stores the tissue images as separate elements of the
    SpatialData object (sdata.images[...]) instead of inside adata.uns["spatial"].
    squidpy's sq.pl.spatial_scatter() needs the classic structure, so we rebuild it
    from the hires/lowres images and from scalefactors_json.json.
    """
    hires = sdata.images[f"{sample}_hires_image"].transpose("y", "x", "c").to_numpy()
    lowres = sdata.images[f"{sample}_lowres_image"].transpose("y", "x", "c").to_numpy()
    with open(f"{data_folder}/{sample}/spatial/scalefactors_json.json") as f:
        scalefactors = json.load(f)
    adata.uns["spatial"] = {
        sample: {
            "images": {"hires": hires, "lowres": lowres},
            "scalefactors": scalefactors,
            "metadata": {},
        }
    }


def plot_spatial_qc(adatas, samples):
    """For each sample (one row): H&E image | n_genes_by_counts | total_counts on the tissue.

    Requires attach_spatial_uns() to have been run on every AnnData in `adatas`.
    """
    import squidpy as sq

    fig, axes = plt.subplots(len(samples), 3, figsize=(15, 5 * len(samples)))
    for row, s in enumerate(samples):
        adata = adatas[s]
        img = adata.uns["spatial"][s]["images"]["hires"]
        axes[row, 0].imshow(img)
        axes[row, 0].set_title(f"{s} - tissue image")
        axes[row, 0].axis("off")
        sq.pl.spatial_scatter(adata, color="n_genes_by_counts", ax=axes[row, 1], title=f"{s} - n_genes_by_counts")
        sq.pl.spatial_scatter(adata, color="total_counts", ax=axes[row, 2], title=f"{s} - total_counts")
    fig.tight_layout()


def relink_table(sdata, adata, samples, region_key="region", instance_key="spot_id"):
    """Re-declare the link between the AnnData table and the Visium spots (shapes).

    A SpatialData table must say which shapes element each row annotates
    (region / region_key / instance_key). This information is stored in
    adata.uns["spatialdata_attrs"] and can be lost when the table is copied or
    modified. This function repairs it, puts the table back in `sdata` and
    returns the repaired AnnData. Always do `adata = eh.relink_table(...)`.
    """
    from spatialdata.models import TableModel

    adata = TableModel.parse(
        adata, region=list(samples), region_key=region_key, instance_key=instance_key
    )
    sdata.tables["table"] = adata
    return adata


def spatial_row(sdata, color, axes, samples, vmax=None, cmap="viridis"):
    """Draw one variable (an .obs column or a gene) on the tissue sections, one per axis.

    sdata : SpatialData object whose "table" contains `color`
    axes  : one matplotlib Axes per sample (e.g. from plt.subplots(1, 4))
    vmax  : if given, all panels share the colour scale [0, vmax] -> comparable panels
    cmap  : colormap (use cmap=None for categorical variables such as clusters)
    """
    import spatialdata_plot  # noqa: F401  (registers the .pl accessor)

    kwargs = {}
    if cmap is not None:
        kwargs["cmap"] = cmap
    if vmax is not None:
        kwargs["norm"] = Normalize(vmin=0, vmax=vmax, clip=True)
    for ax, s in zip(axes, samples):
        sdata.pl.render_images(f"{s}_hires_image").pl.render_shapes(
            s, color=color, table_name="table", method="matplotlib", **kwargs
        ).pl.show(coordinate_systems=f"{s}_downscaled_hires", ax=ax, title=s)


# ----------------------------------------------------------------------------
# Markers and tables
# ----------------------------------------------------------------------------
def top_markers_table(adata, key, n=5):
    """Top-n marker genes per cluster from a rank_genes_groups result, as a tidy table.

    Genes are ranked by test score within each cluster, then displayed by
    decreasing log fold change. Returns a DataFrame (group, names, logfoldchanges).
    """
    top = (
        sc.get.rank_genes_groups_df(adata, group=None, key=key)
        .sort_values("scores", ascending=False)
        .groupby("group", observed=True)
        .head(n)
    )
    top = top.assign(_g=pd.to_numeric(top["group"], errors="coerce"))
    top = top.sort_values(by=["_g", "group", "logfoldchanges"], ascending=[True, True, False])
    return top[["group", "names", "logfoldchanges"]].reset_index(drop=True)


def markers_dict(marker_table, dataset, granularity, adata=None):
    """Turn the paper's marker table into a dict {cell type: [genes]} for sc.pl.dotplot.

    marker_table : DataFrame read from paper_marker_genes.csv
    dataset      : "scRNA-seq" or "ST-seq"
    granularity  : "major" (broad cell types) or the finer level available in the file
    adata        : if given, genes that are not in adata.var_names are dropped
                   (and reported), otherwise the dotplot would fail.
    """
    sel = marker_table[(marker_table["dataset"] == dataset) & (marker_table["granularity"] == granularity)]
    out = {}
    for _, row in sel.iterrows():
        genes = str(row["marker_genes"]).split(";")
        if adata is not None:
            present = [g for g in genes if g in adata.var_names]
            if len(present) < len(genes):
                print(f"{row['cell_type']}: missing genes {sorted(set(genes) - set(present))}")
            genes = present
        if genes:
            out[row["cell_type"]] = genes
    return out


# ----------------------------------------------------------------------------
# Reading the Cell2location results
# ----------------------------------------------------------------------------
def plot_cluster_by_state(q05, adata, cluster_key="leiden", name_key="spot_type"):
    """Mean estimated abundance of each cell state in each ST cluster, as two heatmaps.

    q05 : DataFrame (spots x cell states) of Cell2location q05 abundances
    Left  - composition   : each row (cluster) is normalised to sum to 1
    Right - concentration : each column (cell state) is rescaled to [0, 1]
    Returns (composition, concentration) as DataFrames.
    """
    label = adata.obs[cluster_key].astype(str) + " · " + adata.obs[name_key].astype(str)
    mean_by_cluster = q05.groupby(label.values).mean()
    mean_by_cluster = mean_by_cluster.loc[
        sorted(mean_by_cluster.index, key=lambda x: int(x.split(" · ")[0]))
    ]
    composition = mean_by_cluster.div(mean_by_cluster.sum(axis=1), axis=0)
    concentration = (mean_by_cluster - mean_by_cluster.min()) / (
        mean_by_cluster.max() - mean_by_cluster.min()
    )

    fig, axes = plt.subplots(1, 2, figsize=(20, 5))
    sns.heatmap(composition, cmap="viridis", ax=axes[0], cbar_kws={"label": "fraction of estimated cells"})
    axes[0].set_title("What is each ST cluster made of?")
    sns.heatmap(concentration, cmap="magma", ax=axes[1], cbar_kws={"label": "scaled abundance"})
    axes[1].set_title("Where is each cell state concentrated?")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.tight_layout()
    plt.show()
    return composition, concentration


def compare_state_proportions(q05, adata, adata_sc, samples,
                              sc_to_st=None):
    """Cell-state proportions per sample: Cell2location (ST) next to scRNA-seq.

    q05      : DataFrame (spots x cell states)
    adata    : ST AnnData (needs .obs["sample"])
    adata_sc : scRNA-seq AnnData (needs .obs["sample"] in D0/D1/D7/D30 and .obs["cell_state"])
    sc_to_st : dict translating scRNA-seq timepoints to ST sample names
    Prints the mean estimated number of cells per spot, draws the two stacked bar plots
    and returns (st_prop, sc_prop).
    """
    if sc_to_st is None:
        sc_to_st = {"D0": "Skin", "D1": "Wound1", "D7": "Wound7", "D30": "Wound30"}

    mean_by_sample = q05.groupby(adata.obs["sample"].values, observed=True).mean().loc[samples]
    st_prop = mean_by_sample.div(mean_by_sample.sum(axis=1), axis=0)
    print("Estimated cells per spot (mean of the sum of q05):")
    print(mean_by_sample.sum(axis=1).round(1))

    sc_prop = (
        pd.crosstab(adata_sc.obs["sample"], adata_sc.obs["cell_state"], normalize="index")
        .rename(index=sc_to_st)
        .loc[samples]
        .reindex(columns=st_prop.columns, fill_value=0)
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
    st_prop.plot(kind="bar", stacked=True, colormap="tab20", ax=axes[0], legend=False)
    axes[0].set_title("Spatial data: Cell2location (q05)")
    sc_prop.plot(kind="bar", stacked=True, colormap="tab20", ax=axes[1])
    axes[1].set_title("scRNA-seq: cell states in dissociated cells")
    axes[1].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8)
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("proportion")
        ax.tick_params(axis="x", rotation=0)
    fig.tight_layout()
    plt.show()
    return st_prop, sc_prop
