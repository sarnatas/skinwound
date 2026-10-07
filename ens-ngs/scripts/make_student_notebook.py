#!/usr/bin/env python3
"""
Build Analysis_students.ipynb from the instructor notebook Analysis.ipynb.

    python make_student_notebook.py [teacher.ipynb] [student.ipynb]

The instructor notebook is only READ. Markdown explanations that are kept are copied
verbatim from it (with a sanity check on the expected content, so that if the instructor
notebook is reorganised the script stops instead of silently producing a wrong notebook).
Everything else (task banners, scaffolds, hidden-helper calls) is written here.

Three kinds of cells in the student notebook, always announced by a banner comment:
    RUN        run it as is (often a call to a function hidden in ens_helpers.py)
    READ       tricky code that is shown in full: read it, run it, understand it
    YOUR TURN  the student writes the code (scaffold with ____ blanks or TODO list)
"""
import json
import sys
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEACHER = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "notebooks" / "Analysis.ipynb"
STUDENT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE.parent / "notebooks" / "Analysis_students.ipynb"

nb_t = json.load(open(TEACHER))
T = nb_t["cells"]


def S(i):
    """Source of teacher cell i."""
    return "".join(T[i]["source"])


# ----------------------------------------------------------------------------
# cell builders
# ----------------------------------------------------------------------------
cells = []
_task = {"n": 0}


def _mk(kind, src, tags):
    c = {"cell_type": kind, "id": uuid.uuid4().hex[:8], "metadata": {"tags": tags}, "source": src.strip("\n").splitlines(keepends=True)}
    if kind == "code":
        c["execution_count"] = None
        c["outputs"] = []
    cells.append(c)


def MD(src, tags=None):
    _mk("markdown", src, tags or [])


def CP(i, must=None, replace=None, tags=None):
    """Copy markdown cell i of the teacher notebook (verbatim, optional replace pairs)."""
    s = S(i)
    assert T[i]["cell_type"] == "markdown", f"cell {i} is not markdown"
    if must:
        assert must in s, f"teacher cell {i} changed: expected {must!r}"
    for a, b in (replace or []):
        assert a in s, f"teacher cell {i}: {a!r} not found"
        s = s.replace(a, b)
    MD(s, tags)


def banner(kind, title, lines=()):
    bar = "# " + "=" * 74
    if kind == "YOUR TURN":
        _task["n"] += 1
        head = f"# YOUR TURN  |  Task {_task['n']}  -  {title}"
    else:
        head = f"# {kind}  |  {title}"
    out = [bar, head]
    out += [("# " + l).rstrip() for l in lines]
    out.append(bar)
    return "\n".join(out) + "\n"


def RUN(title, code, lines=()):
    _mk("code", banner("RUN", title, lines) + code, ["run"])


_PLUMBING = ("Save", "Put the filtered", "Reload", "Load your", "Load the Cell2location")


def READ(title, code, lines=()):
    """Tricky code shown in full. Pure plumbing (save/load checkpoints) is labelled RUN instead."""
    if title.startswith(_PLUMBING):
        return RUN(title, code, lines)
    _mk("code", banner("READ", title, lines) + code, ["read"])


def TURN(title, code, lines=()):
    _mk("code", banner("YOUR TURN", title, lines) + code, ["your-turn"])


def Q(text):
    _task.setdefault("q", 0)
    _task["q"] += 1
    MD(f"**Question {_task['q']}.** {text}\n\n*Your answer (double-click to edit):* ", ["question"])


def REVEAL(title, body):
    MD(f"<details>\n<summary><b>{title}</b> (try on your own first, then click)</summary>\n\n{body.strip()}\n\n</details>", ["reveal"])


# ============================================================================
# 0. INTRODUCTION
# ============================================================================
MD("""
# Spatial transcriptomics of human skin wound healing — student notebook

**ENS de Lyon, M2 — NGS practicals.**
Data: Visium spatial transcriptomics (GEO [GSE241124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124)) and scRNA-seq (GEO [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132)) from *Spatiotemporal single-cell roadmap of human skin wound healing* (Liu et al., Cell Stem Cell 2025).

### What you will do
1. Quality control and filtering of Visium spots and genes, normalisation, dimensionality reduction and clustering.
2. Think about **batch effects**, and when *not* to correct them.
3. Analyse the matching scRNA-seq data of the same donor and **annotate cell types** from marker genes.
4. **Annotate the spatial clusters** and **deconvolve** the spots with Cell2location.
5. Read the results: *where* are the cell types, *when*, and *next to whom*?

### How this notebook works
Every code cell starts with a banner that tells you what is expected from you:

| Banner | What it means |
|---|---|
| `RUN` | Just run the cell. Often it calls a function from `ens_helpers.py` that hides repetitive plumbing code; the banner says what the function does (use `help(eh.function_name)` for details). |
| `READ` | Code that is a bit tricky or new. It is shown in full: **run it, but also read it** and be ready to explain what it does. |
| `YOUR TURN` | **You write the code.** Blanks are marked `____`; the banner lists what to do and what you should obtain. |

**Question** cells are for you: write a short answer, we will discuss them together. Boxes titled *click to reveal* contain the discussion — try before you open them!

### Rules of the game
- Run the cells **in order**. If a cell fails, read the error message first.
- Your own results are saved in `my_results/`. The instructor's checkpoints in `data/processed_data/` are a **safety net**: use them only if you are stuck (they are commented out in the notebook). If you load them, remember that **cluster numbers may differ from yours**: redo the annotation dictionaries with *your* clusters.
- Analysis choices matter more than code. Whenever you pick a threshold or a parameter, be ready to justify it with a plot.
""", ["intro"])

# ============================================================================
# 1. SETUP
# ============================================================================
RUN("Libraries", """import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import scanpy as sc
import squidpy as sq
import anndata as ad
import spatialdata as sd
import spatialdata_io
import spatialdata_plot
import seaborn as sns

from scipy.stats import spearmanr
from matplotlib.colors import Normalize

import ens_helpers as eh   # helper functions of this practical (ens_helpers.py, next to this notebook)
""", ["Nothing to modify. If `import ens_helpers` fails, the file ens_helpers.py is not in the same folder as this notebook."])

RUN("Paths and constants  (the only cell where you must edit something)", '''BASE = "/path/to/course/folder"      # <-- EDIT: the folder that contains the `data/` directory

DATA        = f"{BASE}/data/GSE241124_RAW"                 # Visium: one folder per sample
SC_DATA     = f"{BASE}/data/GSE241132_RAW"                 # scRNA-seq: one folder per sample
PROC        = f"{BASE}/data/processed_data"                # instructor checkpoints (read only, safety net)
MARKERS_CSV = f"{BASE}/data/paper_marker_genes.csv"        # marker genes reported in the paper
PAPER_META  = f"{SC_DATA}/GSE241132_cell_metadata.txt"     # per-cell metadata published by the authors
Q05_FILE    = f"{PROC}/q05_cell_abundance.csv"             # Cell2location output (computed by the instructors)

MY = f"{BASE}/my_results"                                  # YOUR checkpoints are written here
os.makedirs(MY, exist_ok=True)

SAMPLES = ["Skin", "Wound1", "Wound7", "Wound30"]          # the 4 Visium sections (same donor, 4 timepoints)
METRICS = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
TITLES = {
    "total_counts": "Total UMIs per spot",
    "n_genes_by_counts": "Detected genes per spot",
    "pct_counts_mt": "% mitochondrial counts per spot",
}

sns.set_style("whitegrid")
''', ["Set BASE to the course folder. Everything else is derived from it.",
      "If a path error appears later on, this is the first place to check."])

# ============================================================================
# 2. SPATIAL TRANSCRIPTOMICS
# ============================================================================
CP(2, "Spatial Transcriptomics data analysis")
MD("""
### Loading the Visium data

We analyse the sections of **one donor** at four timepoints: healthy skin (`Skin`) and wound-edge biopsies at day 1, 7 and 30 (`Wound1`, `Wound7`, `Wound30`). Each section is read into a [SpatialData](https://spatialdata.scverse.org) object, which keeps together the expression table, the spot coordinates (shapes) and the tissue images.
""")

READ("Load one SpatialData object per sample", '''sdatas = {
    s: spatialdata_io.visium(f"{DATA}/{s}", dataset_id=s)
    for s in SAMPLES
}
''', ["Nothing to modify. `spatialdata_io.visium` reads a Space Ranger output folder."])

TURN("Inspect a SpatialData object", """# What to do:
#   1. Display one object, e.g. sdatas["Skin"]. Which elements does it contain (images, shapes, tables)?
#   2. Look at its table: how many spots and genes? What is in .obs? Where are the spot coordinates stored (.obsm)?
# Expected: you can answer the question below.

# write your code here
""")
Q("What does each element of the SpatialData object contain, and what links the expression table to the spots drawn on the image?")

MD("""
### Quality control

We compute, for every spot, three standard QC metrics: the **total number of UMIs**, the **number of detected genes** and the **percentage of mitochondrial counts**. The function `eh.qc_metrics` flags the mitochondrial genes (symbol starting with `MT-`) and calls `sc.pp.calculate_qc_metrics`.
""")

RUN("QC metrics for every sample", '''adatas = {}
for s in SAMPLES:
    table = sdatas[s]["table"].copy()
    eh.qc_metrics(table)            # adds total_counts, n_genes_by_counts, pct_counts_mt to table.obs
    table.obs["sample"] = s
    adatas[s] = table
    print(f"{s}: {table.n_obs} spots, {table.n_vars} genes")
''', ["Expected: one line per sample with the number of spots and genes."])

RUN("QC violin plots, one figure per sample", """eh.plot_qc_violins(adatas, METRICS, TITLES, unit="spots")
""", ["Expected: 4 figures with 3 panels each (UMIs, genes, % mitochondrial).",
      "Look carefully: you will need these plots to choose your filtering thresholds below."])

READ("Combine the four samples into one object and compare them", '''# sd.concatenate merges the SpatialData objects. concatenate_tables=True also concatenates
# the underlying tables, so all spots of all samples end up in one shared table.
sdata_all = sd.concatenate(list(sdatas.values()), concatenate_tables=True)
table_all = sdata_all["table"].copy()      # the merged expression table
eh.qc_metrics(table_all)                    # QC metrics on the merged table
table_all.obs.head()
''', ["Run it and look at the columns of table_all.obs.",
      "Note the `region` column: it is created by spatialdata_io and contains the sample name."])

RUN("QC of the four samples side by side", """eh.plot_qc_combined(table_all.obs, "region", SAMPLES, METRICS, TITLES, "QC combined across samples")
""", ["Expected: 3 panels (one per metric), 4 violins per panel."])
Q("Compare the four samples. Do they have the same sequencing depth? Which sample looks different, and what could explain it?")

CP(15, "spatial_scatter", replace=[("sq.pl.spatial_scatter(adata color=", "sq.pl.spatial_scatter(adata, color=")])
MD("""*In this practical the function that rebuilds this structure (`eh.attach_spatial_uns`) is provided: you do not have to write it.*""")

RUN("Visualise the QC metrics on the tissue", """for s in SAMPLES:
    eh.attach_spatial_uns(adatas[s], sdatas[s], s, DATA)    # rebuild adata.uns["spatial"] for squidpy

eh.plot_spatial_qc(adatas, SAMPLES)
""", ["Expected: a grid with one row per sample: tissue image | detected genes | total UMIs.",
      "Is low quality random in space, or does it follow tissue structures?"])

MD("""
### Filtering spots

Spots with very few detected genes are mostly empty or poorly captured spots. Where to cut is a judgement call that **depends on the sample** (the sequencing depth is not the same): use the per-sample violin plots above.
""")

TURN("Choose the minimum number of detected genes per sample", """# What to do:
#   Replace each None with the minimum number of detected genes (n_genes_by_counts) that a spot must have
#   in that sample, by looking at where the low-quality tail starts in the violin plots.
#   Write one sentence below justifying your choice for each sample.
# Expected: a few % of spots removed per sample, not 50%. Check the printout.

MIN_GENES_PER_SAMPLE = {
    "Skin": None,
    "Wound1": None,
    "Wound7": None,
    "Wound30": None,
}
assert all(v is not None for v in MIN_GENES_PER_SAMPLE.values()), "Fill in a threshold for every sample"

for s in SAMPLES:
    n_before = adatas[s].n_obs
    sc.pp.filter_cells(adatas[s], min_genes=MIN_GENES_PER_SAMPLE[s])
    n_after = adatas[s].n_obs
    print(
        f"{s}: min_genes={MIN_GENES_PER_SAMPLE[s]} -> "
        f"{n_before} to {n_after} spots ({n_before - n_after} removed, "
        f"{(n_before - n_after) / n_before * 100:.1f}%)"
    )
""")
MD("*Your justification (double-click to edit):* ")

MD("""
### Filtering genes, collectively across samples

We now put the four filtered samples into one table and remove genes that are detected in almost no spot.
""")

READ("Concatenate the filtered samples", '''table_all = ad.concat(
    [adatas[s] for s in SAMPLES],
    join="inner",        # only keep the genes present in all samples
    index_unique="-",    # make the spot barcodes unique by appending the sample key
    keys=SAMPLES,        # one key per concatenated object
)
''', ["Nothing to modify. Read the comments: what would `join=\"outer\"` do instead?"])

TURN("Remove genes detected in very few spots", """# What to do:
#   Keep only the genes detected in at least 3 spots (sc.pp.filter_genes, argument min_cells).
# Expected: a printout with the number of genes before/after, the table shape and the spots per sample.

n_genes_before = table_all.n_vars
____                                    # <- your filtering call
print(f"Genes: {n_genes_before} -> {table_all.n_vars} ({n_genes_before - table_all.n_vars} removed)")
print(f"Combined table: {table_all.shape[0]} spots x {table_all.shape[1]} genes")
print(table_all.obs["sample"].value_counts())
""")

MD("""
### Do we need to filter spots on the percentage of mitochondrial counts?

In dissociated single-cell data, a high mitochondrial fraction is a classic sign of dying cells, and cells above ~20% are usually removed. Does the same logic apply to a **tissue section**?
""")
RUN("How many spots would a 20% mitochondrial threshold remove?", """# demo only - NOT applied to the data
n_would_remove = (table_all.obs["pct_counts_mt"] > 20).sum()
print(
    f"Spots with pct_counts_mt > 20%: {n_would_remove} / {table_all.n_obs} "
    f"({n_would_remove / table_all.n_obs * 100:.2f}%)"
)
""", ["Nothing to modify. We only look at the number: no spot is removed."])
Q("Based on this number and on what a Visium spot is, would you filter on `pct_counts_mt` in this dataset? Why (not)? And how could a very active tissue region affect this metric?")
s27 = S(27)
assert "In this dataset that threshold would remove almost nothing" in s27
REVEAL("Discussion: mitochondrial counts in tissue, and the strategy of the paper", s27)

TURN("Remove MALAT1, mitochondrial and hemoglobin genes", """# What to do:
#   Following the paper, remove from the expression table (collectively, not per sample):
#     - MALAT1
#     - all mitochondrial genes (gene symbol starting with "MT-")
#     - the hemoglobin genes listed below (only those that exist in the table)
#   Hint: table_all.var_names.str.startswith(...), .isin(...), and don't forget .copy().
#   PITFALL: do NOT select hemoglobin genes with a prefix like "HB": you would also remove HBEGF,
#   a real wound-healing growth factor (and HBP1, HBS1L). Use the explicit list.
# Expected: a few dozen genes removed at most.

HEMOGLOBIN_GENES = [
    "HBA1", "HBA2", "HBB", "HBD", "HBE1",
    "HBG1", "HBG2", "HBM", "HBQ1", "HBZ",
]

mt_genes = ____            # list of mitochondrial genes present in the table
hb_genes = ____            # hemoglobin genes present in the table
malat1_gene = ____         # ["MALAT1"] if present
genes_to_remove = set(malat1_gene) | set(mt_genes) | set(hb_genes)

print(f"Mitochondrial genes ({len(mt_genes)}): {mt_genes}")
print(f"Hemoglobin genes ({len(hb_genes)}): {hb_genes}")
print(f"Total genes removed: {len(genes_to_remove)}")

n_genes_before = table_all.n_vars
table_all = ____           # your filtered table (a copy)
print(f"Genes: {n_genes_before} -> {table_all.n_vars}")
print(f"Final table_all: {table_all.shape[0]} spots x {table_all.shape[1]} genes")
""")

CP(29, "Put back all the info")
READ("Put the filtered table back in the SpatialData object and save it", '''sdata_all.tables["table"] = table_all      # replaces the original table by the filtered one
sdata_all.write(f"{MY}/spatial_data_filtered.zarr", overwrite=True)
''', ["Nothing to modify. It writes a checkpoint in your my_results/ folder."])

READ("Reload the checkpoint", '''sdata_all = sd.read_zarr(f"{MY}/spatial_data_filtered.zarr")
# --- SAFETY NET: if you are stuck above, load the instructor's checkpoint instead (uncomment) ---
# sdata_all = sd.read_zarr(f"{PROC}/spatial_data_filtered.zarr")

adata = sdata_all["table"]
# spatialdata_io calls the sample column "region", which is misleading here (it is always the same
# wound, at different timepoints). We copy it to a column called "sample", used from now on.
adata.obs["sample"] = adata.obs["region"].astype(str)
print(adata.shape, adata.obs["sample"].unique().tolist())
sdata_all
''', ["Expected: the shape of the table and the 4 sample names, then the SpatialData summary."])

MD("""
### Normalisation, highly variable genes, PCA
""")
READ("Keep the raw counts", '''# keep raw counts: needed later by Cell2location, which requires raw (not normalised) counts
adata.layers["counts"] = adata.X.copy()
''', ["Nothing to modify. Never forget this step before normalising!"])
CP(37, "We normalize with a standard")

TURN("Normalisation, highly variable genes, PCA", """# What to do (3 steps, in this order):
#   1. Normalise each spot to 10,000 total counts, then log-transform      -> sc.pp.normalize_total, sc.pp.log1p
#   2. Select 2,000 highly variable genes, batch-aware (batch_key="sample") -> sc.pp.highly_variable_genes
#   3. PCA with 30 components (it uses the highly variable genes by default) -> sc.tl.pca
# Expected: adata.var["highly_variable"] has ~2000 True values, adata.obsm["X_pca"] has 30 columns.

____
____
____

print(adata.obsm["X_pca"].shape, int(adata.var["highly_variable"].sum()))
""")

MD("""
### Uncorrected UMAP, colored by sample
""")
READ("Neighbors graph and UMAP, without batch correction", '''sc.pp.neighbors(adata, n_pcs=30, key_added="neighbors_uncorrected")
sc.tl.umap(adata, neighbors_key="neighbors_uncorrected")
# Each new UMAP overwrites adata.obsm["X_umap"]: we store a copy of this one under another name.
adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()

sc.pl.embedding(adata, basis="X_umap_uncorrected", color="sample", title="uncorrected")
''', ["Nothing to modify. Note `key_added`: it lets us keep several neighbor graphs in the same object."])

MD("""
### Batch effects: to correct or not to correct?

The four samples separate on the UMAP above. Four different datasets put together, a clear separation… the textbook reflex is to *remove the batch effect*.
""")
Q("Before correcting anything: what exactly does the variable `sample` represent in this dataset? Is it a technical batch, a biological condition, or both? What would you lose by correcting for it?")
s43 = S(43)
split_marker = "Now, just to show you what are 3 among the most common tools"
assert split_marker in s43
part1, part2 = s43.split(split_marker)
REVEAL("Discussion: batch and biology are confounded here", part1.replace("### Note about batch effect removal", "").strip())
MD(split_marker + part2)

READ("BBKNN: correct the neighbor graph", S(44), ["Nothing to modify. Run it and look at the UMAP."])
READ("Scanorama: correct the expression matrix", S(45), ["Nothing to modify. It is the slowest of the three: be patient (and read the comments)."])
READ("Harmony: correct the PCA embedding", S(46), ["Nothing to modify. Note the workaround for the broken scanpy wrapper: harmonypy is called directly."])
Q("The three UMAPs are very different. What does each method correct (graph, embedding, expression)? Do you trust one more than the others? In which situation would correcting for batch be the right thing to do here?")

CP(47, "Leiden clustering on uncorrected data")
TURN("Leiden clustering on the uncorrected data", """# What to do:
#   1. Run the Leiden algorithm (sc.tl.leiden) on the UNCORRECTED neighbor graph (the one stored as
#      "neighbors_uncorrected"), resolution 0.5, and store the result in adata.obs["leiden"].
#   2. Print the number of spots per cluster.
#   3. Plot the clusters on the uncorrected UMAP (basis "X_umap_uncorrected").
# Expected: about 7 clusters.

____
print(adata.obs["leiden"].value_counts())

____
""")

READ("Save a checkpoint with the clusters", '''sdata_all.tables["table"] = adata
sdata_all.write(f"{MY}/spatial_data_clustered.zarr", overwrite=True)
''', ["Nothing to modify. Next time, you can restart from here."])

MD("""
### Clusters on the tissue

To plot on the tissue we have to make sure that the table is correctly linked to the spots (shapes) of the SpatialData object. This bookkeeping is done by `eh.relink_table` (it re-declares the table annotation, which can be lost when the table is modified).
""")
RUN("Plot the Leiden clusters on the tissue sections", """adata = eh.relink_table(sdata_all, adata, SAMPLES)

fig, axes = plt.subplots(1, len(SAMPLES), figsize=(5 * len(SAMPLES), 5))
eh.spatial_row(sdata_all, "leiden", axes, SAMPLES, cmap=None)
fig.tight_layout()
plt.show()
""", ["Expected: one panel per sample, spots coloured by cluster.",
      "Do the clusters follow the tissue structures you see in the H&E image?"])
Q("Describe how the clusters are organised in space. Do some clusters appear only in certain samples? What does it tell you?")

CP(53, "Locating the wound within each section")
RUN("Wound-edge marker genes on the tissue", """wound_edge_markers = ["KRT6A", "KRT16", "KRT17", "S100A8"]

for gene in wound_edge_markers:
    fig, axes = plt.subplots(1, len(SAMPLES), figsize=(5 * len(SAMPLES), 5))
    eh.spatial_row(sdata_all, gene, axes, SAMPLES)
    fig.suptitle(gene)
    fig.tight_layout()
    plt.show()
""", ["Expected: 4 figures (one per gene), 4 panels each."])
Q("Where is the wound edge in each section? Do the four genes agree with each other? Which timepoint shows the strongest signal, and is it surprising?")

# ============================================================================
# 3. SINGLE CELL
# ============================================================================
CP(59, "Single-cell data analysis")
CP(60, "Why we cluster scRNA-seq and ST-seq independently")

MD("""
### Loading the scRNA-seq data of the same donor

We use the scRNA-seq samples of donor **PWH27** (the same donor as the Visium sections), at four timepoints: D0 (healthy skin) and D1, D7, D30 after wounding.
""")
READ("Load and concatenate the four timepoints", '''donor = "PWH27"
timepoints = ["D0", "D1", "D7", "D30"]

adatas_sc = {}
for tp in timepoints:
    a = sc.read_10x_mtx(f"{SC_DATA}/{donor}{tp}")
    a.obs["donor"] = donor
    a.obs["sample"] = tp
    a.obs_names = [f"{bc}-{tp}" for bc in a.obs_names]      # make barcodes unique across timepoints
    adatas_sc[tp] = a

adata_sc = ad.concat(list(adatas_sc.values()), join="inner")
adata_sc
''', ["Nothing to modify. Note how the barcodes are renamed: we will need this convention later,",
      "when we match our cells with the annotation published by the authors."])

RUN("QC metrics and violin plots (same as for the spatial data)", '''eh.qc_metrics(adata_sc)

SC_SAMPLES = ["D0", "D1", "D7", "D30"]
SC_TITLES = {
    "total_counts": "Total UMIs per cell",
    "n_genes_by_counts": "Detected genes per cell",
    "pct_counts_mt": "% mitochondrial counts per cell",
}

eh.plot_qc_violins({s: adata_sc[adata_sc.obs["sample"] == s] for s in SC_SAMPLES}, METRICS, SC_TITLES, unit="cells")
eh.plot_qc_combined(adata_sc.obs, "sample", SC_SAMPLES, METRICS, SC_TITLES, "QC combined across samples (donor PWH27)")
''', ["Expected: 4 figures (one per timepoint) and a combined figure."])

TURN("Filter low-quality cells", """# What to do: apply the same cell-level filters as the paper
#   - at least 500 detected genes per cell          (sc.pp.filter_cells, min_genes)
#   - at least 1000 UMIs per cell                   (sc.pp.filter_cells, min_counts)
#   - less than 20% mitochondrial counts            (boolean mask on adata_sc.obs["pct_counts_mt"], then .copy())
# Note: here (dissociated cells), the mitochondrial filter IS appropriate. Why is it different from the tissue?
# Expected: the number of cells decreases at each step. Print it.

print("before:", adata_sc.n_obs)
____
____
____
print("after:", adata_sc.n_obs)
""")

MD("""
### Doublets

Two cells captured in the same droplet give a mixed profile that can look like a new cell type. **Scrublet** simulates doublets from the data and scores each cell by its similarity to them.
""")
READ("Detect and remove doublets with Scrublet", '''sc.pp.scrublet(adata_sc, batch_key="sample")       # run per sample: doublets only form within one droplet batch
print(adata_sc.obs["predicted_doublet"].value_counts())

n_before = adata_sc.n_obs
adata_sc = adata_sc[~adata_sc.obs["predicted_doublet"]].copy()
print(f"Doublet removal: {n_before} -> {adata_sc.n_obs} cells")
''', ["Nothing to modify. Why do we run it per sample (batch_key)?"])

TURN("Gene filtering", """# What to do: same strategy as for the spatial data
#   1. Remove MALAT1, the mitochondrial genes and the hemoglobin genes (reuse HEMOGLOBIN_GENES).
#   2. Then keep only the genes detected in at least 10 cells (sc.pp.filter_genes, min_cells).
# Expected: a few dozen genes removed in step 1, a few thousand in step 2.

n_genes_before = adata_sc.n_vars
____
____
print(f"Genes: {n_genes_before} -> {adata_sc.n_vars}")
""")

TURN("Normalisation, PCA, UMAP (uncorrected)", """# What to do: repeat the pipeline you wrote for the spatial data, on adata_sc:
#   1. normalise to 10,000 counts + log1p
#   2. 2,000 highly variable genes, batch-aware (batch_key="sample")
#   3. PCA (30 components)
#   4. neighbors graph (30 PCs), then UMAP
#   5. plot the UMAP coloured by "sample"
# Expected: a UMAP in which the timepoints are partially separated.

adata_sc.layers["counts"] = adata_sc.X.copy()      # raw counts, needed later by Cell2location (already written for you)

____
""")

READ("Save a checkpoint, and a batch correction demo", '''adata_sc.write_h5ad(f"{MY}/sc_data_PWH27_processed.h5ad")

# Batch correction demo with Harmony (same logic and same caveats as for the spatial data)
ho_sc = harmonypy.run_harmony(adata_sc.obsm["X_pca"], adata_sc.obs, ["sample"])
adata_sc.obsm["X_pca_harmony"] = ho_sc.Z_corr
sc.pp.neighbors(adata_sc, use_rep="X_pca_harmony", key_added="harmony")
sc.tl.umap(adata_sc, neighbors_key="harmony")
sc.pl.umap(adata_sc, neighbors_key="harmony", color="sample", title="harmony (single donor)")

# The Harmony UMAP just overwrote adata_sc.obsm["X_umap"]. We continue with the UNCORRECTED data,
# as decided: recompute the UMAP on the default ("neighbors") graph.
sc.tl.umap(adata_sc)
sc.pl.umap(adata_sc, color="sample", title="uncorrected (recomputed)")
''', ["Nothing to modify. Same reasoning as before: with a single donor, timepoint and batch are confounded."])

MD("""
### Choosing the clustering resolution

The resolution of Leiden clustering controls how many clusters we get. Instead of guessing, **RESOLUTE** scans a range of resolutions and scores each of them with several criteria (Calinski–Harabasz, BIC, and bootstrap stability).
""")
READ("RESOLUTE with the Calinski-Harabasz score (and bootstrap stability)", S(76), ["Nothing to modify. This can take a few minutes. Look at the plots it produces."])
READ("RESOLUTE with the BIC score", S(77), ["Nothing to modify. Look at the BIC curve: does it show a clear optimum?"])
Q("Compare the three criteria (Calinski-Harabasz, BIC, stability). Do they agree? Which resolution would you choose, and why? Remember that the resolution is a working choice, validated afterwards by the biology.")
s78 = S(78)
assert "resolution 0.2" in s78
REVEAL("What we chose in our run", s78.replace("## Choosing the clustering resolution", "").strip())

TURN("Final Leiden clustering", """# What to do: set CHOSEN_RESOLUTION to the value you decided, then run the cell.
CHOSEN_RESOLUTION = None
assert CHOSEN_RESOLUTION is not None, "Choose a resolution first"

sc.tl.leiden(
    adata_sc,
    neighbors_key="neighbors",
    resolution=CHOSEN_RESOLUTION,
    key_added="leiden",
    flavor="igraph",
    directed=False,
)
print(adata_sc.obs["leiden"].value_counts())

fig, axes = plt.subplots(1, 2, figsize=(15, 5))
sc.pl.embedding(adata_sc, basis="X_umap", color="sample", ax=axes[0], show=False, title="colored by sample")
sc.pl.embedding(adata_sc, basis="X_umap", color="leiden", ax=axes[1], show=False, title=f"Leiden (resolution={CHOSEN_RESOLUTION})")
fig.tight_layout()
plt.show()
""", ["Nothing else to write: just fill CHOSEN_RESOLUTION. Expected: ~15-20 clusters at the resolution we used."])

MD("""
### Marker genes of each cluster
""")
TURN("Differential expression between clusters", """# What to do:
#   1. One-vs-rest differential expression between the Leiden clusters of adata_sc, Wilcoxon test,
#      results stored under key "rank_genes_leiden"                          -> sc.tl.rank_genes_groups
#   2. Compute the dendrogram of the clusters on the PCA                     -> sc.tl.dendrogram
#   3. Dotplot of the top 5 genes of each cluster, scaled per gene            -> sc.pl.rank_genes_groups_dotplot
#      (arguments: groupby, key, n_genes, standard_scale="var")
# Expected: a dotplot with ~5 genes per cluster.

____
____
____
""")

RUN("The marker genes published in the paper", """marker_table = pd.read_csv(MARKERS_CSV)
display(marker_table)

# dict {major cell type: [marker genes]} for the scRNA-seq data (genes missing from our data are dropped)
sc_major_markers = eh.markers_dict(marker_table, "scRNA-seq", "major", adata_sc)

sc.pl.dotplot(adata_sc, sc_major_markers, groupby="leiden", standard_scale="var", dendrogram=True)
""", ["Expected: the table of the paper's markers, then a dotplot (cell types of the paper x our clusters).",
      "Which clusters light up for which cell type?"])

READ("Save a checkpoint with the clusters", '''adata_sc.write_h5ad(f"{MY}/adata_sc_leiden.h5ad")
''', ["Nothing to modify."])

CP(85, "Cell type assignment - from clusters to cell types",
   replace=[("(here, 17 clusters)", "(about 15-20 clusters with the resolution we used)")])

RUN("Top 5 marker genes per cluster, as a readable table", """top5 = eh.top_markers_table(adata_sc, key="rank_genes_leiden", n=5)
print(top5.to_string())
""", ["Expected: 5 genes per cluster, ranked by log fold change.", "Keep this table at hand for the annotation."])

MD("""
### Your annotation

For **each** Leiden cluster you now have what you need to decide which cell type it is:
the top marker genes (table above), the dotplot of the paper's markers (above), and your knowledge of skin biology.
A reasonable workflow:

1. Read the top markers of the cluster. Do they match a known cell type? (Be careful: genes such as *RPS/RPL*, *B2M*, *TMSB4X*, or metallothioneins *MT1X/MT1G* are not specific; a cluster dominated by them may be low quality, stressed, or simply unspecific.)
2. Look at the paper's marker dotplot for the same cluster.
3. For ambiguous clusters, run a **targeted dotplot** with the markers of the candidate cell types (see the next task).
4. Write your label. It is perfectly normal that **several clusters map to the same cell type** (different states of the same population); and rare populations (for instance Schwann cells) may not form a cluster of their own at all.

The finer descriptions (basal, differentiated, wound-activated...) are *interpretations based on your knowledge* of skin biology: they are working hypotheses, not validated annotations.
""")

TURN("Targeted check of ambiguous clusters", """# What to do (optional but recommended for clusters you are unsure about):
#   build a small dict {"group name": [marker genes]} with the genes of the candidate cell types
#   (e.g. T/NK cells, or any population you suspect), and plot it per cluster with sc.pl.dotplot.
# Expected: a dotplot in which one cluster stands out for the genes you chose.

targeted_markers = {
    # "T/NK markers": ["CD3D", "CD3E", "NKG7"],     # <- example of the format
}

____
""")

TURN("Assign a cell type to every cluster", '''# What to do: complete the dictionary {cluster number (as a string): cell type} for ALL clusters.
# Suggested labels (the major cell types of the paper that are expected in this subset):
#   "Keratinocytes", "Fibroblasts", "Myeloid cells", "Lymphoid cells", "Mast cells",
#   "Melanocytes", "Pericytes and smooth muscle cells", "Endothelial cells"
# You may use another label if you can justify it.
# Expected: a UMAP coloured by cell type, with a handful of well-separated groups.

cluster_to_celltype = {
    # "0": "...",
}

missing = sorted(set(adata_sc.obs["leiden"].unique()) - set(cluster_to_celltype))
assert not missing, f"Clusters without annotation: {missing}"

adata_sc.obs["cell_type"] = adata_sc.obs["leiden"].map(cluster_to_celltype).astype("category")
sc.pl.embedding(adata_sc, basis="X_umap", color="cell_type", title="Cell type annotation")
''')

CP(92, "Validating the annotation against the published dataset")

READ("Read the metadata published by the authors", '''paper_meta = pd.read_csv(PAPER_META, sep="\\t", index_col=0)

print("Paper metadata barcode format (first 5):")
print(paper_meta.index[:5].tolist())
print("\\nYour adata_sc barcode format (first 5):")
print(adata_sc.obs_names[:5].tolist())
print("\\nPaper metadata columns:", paper_meta.columns.tolist())
''', ["Nothing to modify. Compare the two barcode formats: they differ, hence the next cell."])

READ("Match our cells to the published annotation and cross-tabulate", S(93),
     ["Nothing to modify. Read the function to_paper_barcode.", "The cross-tabulation: rows = YOUR clusters, columns = the cell types assigned by the authors."])

TURN("Do you confirm or correct your annotation?", """# What to do:
#   Look at the cross-table above, cluster by cluster. For each cluster, which cell type did the authors
#   assign to most of its cells? Does it agree with your label? (Names differ slightly: "Keratinocyte"
#   vs "Keratinocytes", etc.)
#   Optionally, quantify the agreement with ONE number of your choice (e.g. the fraction of cells in
#   each cluster that falls in its majority published type), using the cross-table.
#   If you find a clear disagreement, go back to your dictionary above, correct it, and re-run.

# write your code here
""")
Q("Which clusters, if any, disagree with the published annotation? Was it a mistake in your annotation, or is the published label a different level of granularity?")

READ("Save a checkpoint with the annotation", '''adata_sc.write_h5ad(f"{MY}/adata_sc_leiden_annotated.h5ad")

# --- SAFETY NET: instructor's annotated object (uncomment if you are stuck; cluster numbers are the instructor's!) ---
# adata_sc = sc.read_h5ad(f"{PROC}/adata_sc_leiden_res0.2_annotated.h5ad")
''', ["Nothing to modify."])

MD("""
### A finer label for the deconvolution: cell states

Cell2location needs a reference with **one label per cell state**, finer than the major cell type (for instance, *basal*, *differentiated* or *wound-activated* keratinocytes). The next step labels each cluster with a **state**.

**Important:** the deconvolution that we will read later was computed by the instructors with a fixed list of 17 cell states. Your labels must come from this list, so that your results can be compared with it.
""")
TURN("Assign a cell state to every cluster", '''# What to do: complete {cluster number (string): cell state}, using ONLY the labels of this list.
# Think about which marker genes tell the states apart (you may need extra targeted dotplots).
CELL_STATES = [
    "Keratinocytes_basal_A", "Keratinocytes_basal_B", "Keratinocytes_differentiated",
    "Keratinocytes_suprabasal", "Keratinocytes_ribosomal_stress", "Keratinocytes_wound_activated",
    "Myeloid_MHCII_high", "Macrophages", "Lymphoid_T_NK", "Mast_cells",
    "Fibroblasts_CFD", "Fibroblasts_myofibroblast_like", "Fibroblasts_ECM_high",
    "Melanocytes", "Pericytes_SMC", "Endothelial_lymphatic", "Endothelial_vascular",
]
# Several clusters may share a state; some states may be absent from your clustering.
# Hint: the state must be coherent with the cell type you gave to the same cluster.

cluster_to_state = {
    # "0": "...",
}

missing = sorted(set(adata_sc.obs["leiden"].unique()) - set(cluster_to_state))
unknown = sorted(set(cluster_to_state.values()) - set(CELL_STATES))
assert not missing, f"Clusters without a state: {missing}"
assert not unknown, f"Labels not in CELL_STATES: {unknown}"

adata_sc.obs["cell_state"] = adata_sc.obs["leiden"].map(cluster_to_state).astype("category")
print(adata_sc.obs["cell_state"].value_counts())
''')

# ============================================================================
# 4. BACK TO THE SPATIAL DATA: ANNOTATE SPATIAL CLUSTERS
# ============================================================================
CP(99, "Coming back to the Visium 10x data")

READ("Load your clustered spatial data", '''sdata_all = sd.read_zarr(f"{MY}/spatial_data_clustered.zarr")
# --- SAFETY NET (uncomment if you are stuck; cluster numbers are the instructor's!) ---
# sdata_all = sd.read_zarr(f"{PROC}/spatial_data_clustered.zarr")
adata = sdata_all["table"]
''', ["Nothing to modify."])

TURN("Differential expression between spatial clusters", """# What to do: the same analysis you did for the single-cell clusters, now on the spatial clusters.
#   1. Print how many spots are in each cluster.
#   2. One-vs-rest Wilcoxon DE between the Leiden clusters of `adata`, stored under key "rank_genes_leiden".
#   3. Print the top 5 genes per cluster as a readable table (helper: eh.top_markers_table).
#   4. Dendrogram + dotplot of the top 5 genes per cluster (same functions as for the single cells).
# Expected: a table with 5 genes per cluster and a dotplot.

adata = sdata_all["table"]
____
____
____
____
""")

RUN("The marker genes of the paper on the spatial clusters", """marker_table = pd.read_csv(MARKERS_CSV)
st_markers = eh.markers_dict(marker_table, "ST-seq", "major", adata)

sc.pl.dotplot(adata, st_markers, groupby="leiden", standard_scale="var", dendrogram=True)
""", ["Expected: a dotplot (the paper's spatial marker categories x your spatial clusters)."])

MD("""
#### Spot-type labels: an interpretation

You will now give each spatial cluster a name, from its top DE markers, the paper's ST marker categories, the tissue image and the position of the clusters on the sections (use the maps above). Remember that a spot is a **mixture** of cells: a cluster is *dominated by* a cell type or a tissue structure, it is not a pure population.

Your `spot_type` labels are **working interpretations, not validated annotations**. In particular, think about which clusters you are unsure about: they are good questions for the deconvolution!
""")
TURN("Name the spatial clusters", '''# What to do: complete {cluster number (string): short descriptive name} for ALL spatial clusters.
# Examples of the kind of names you can use (they are NOT the answer, they only show the format):
#   "Suprabasal_epidermis", "Dermis_fibroblasts", "Inflamed_tissue"
# Think about: epidermis vs dermis, wound edge vs intact tissue, immune infiltrate, vessels...

cluster_to_spot_type = {
    # "0": "...",
}

missing = sorted(set(adata.obs["leiden"].unique()) - set(cluster_to_spot_type))
assert not missing, f"Clusters without a name: {missing}"

adata.obs["spot_type"] = adata.obs["leiden"].map(cluster_to_spot_type).astype("category")
adata.obs[["leiden", "spot_type"]].drop_duplicates().sort_values("leiden")
''')
Q("For which clusters are you least sure of your label, and why? What additional information would help you decide?")

# ============================================================================
# 5. CELL2LOCATION
# ============================================================================
s106 = S(106)
cut = s106.index("Last but not least")
MD(s106[:cut].rstrip() + """

**In this practical you do not train the model**: training takes hours on a CPU (a GPU is strongly recommended) and it is better run from a script than from a notebook. The instructors ran it for you and provide the output table, `q05_cell_abundance.csv`. Below, the code that was run is **shown in full** (not executable on purpose): read it and make sure you understand each step. The curious can find the standalone script `run_cell2location.py` in the course repository.
""")

MD("""
#### The code that was run (read only)

**Setup**
```python
import numpy as np
import pandas as pd
import torch
import scvi
import matplotlib.pyplot as plt
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel, Cell2location

N_THREADS = 8                  # number of cores you can use
REF_EPOCHS = 250
VIS_EPOCHS = 5000              # use ~200 for a quick demo

torch.set_num_threads(N_THREADS)
scvi.settings.seed = 0
ACC = "gpu" if torch.cuda.is_available() else "cpu"
```

**Stage 1 - reference signatures from the scRNA-seq (raw counts, labels = `cell_state`)**
```python
adata_ref = adata_sc.copy()
adata_ref.X = adata_ref.layers["counts"].copy()

# Keep informative genes only
selected = filter_genes(adata_ref, cell_count_cutoff=5,
                        cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

RegressionModel.setup_anndata(adata_ref, batch_key="sample", labels_key="cell_state")
mod_ref = RegressionModel(adata_ref)
mod_ref.train(max_epochs=REF_EPOCHS, accelerator=ACC)
adata_ref = mod_ref.export_posterior(
    adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": ACC})

# Average expression signature of each cell state
states = adata_ref.uns["mod"]["factor_names"]
inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{s}" for s in states]].copy()
inf_aver.columns = states
```

**Stage 2 - mapping onto the Visium spots (raw counts of the spatial data)**
```python
adata_vis = adata.copy()
adata_vis.X = adata_vis.layers["counts"].copy()

# Use the genes shared by spatial data and signatures
shared = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, shared].copy()
inf_aver = inf_aver.loc[shared, :].copy()

Cell2location.setup_anndata(adata_vis, batch_key="sample")
mod_vis = Cell2location(adata_vis, cell_state_df=inf_aver,
                        N_cells_per_location=20, detection_alpha=20)
mod_vis.train(max_epochs=VIS_EPOCHS, batch_size=None, train_size=1, accelerator=ACC)

adata_vis = mod_vis.export_posterior(
    adata_vis, sample_kwargs={"num_samples": 1000,
                              "batch_size": mod_vis.adata.n_obs, "accelerator": ACC})
q05 = adata_vis.obsm["q05_cell_abundance_w_sf"].copy()      # the table you will receive
```
""")
Q("Why does Cell2location need *raw counts* (not normalised) for both the reference and the spatial data? What do `N_cells_per_location` and `detection_alpha` represent, and where does the value 20 come from?")
Q("Stage 1 learns a *signature* per cell state, stage 2 only *uses* those signatures. What does this imply for a structure that is present in the tissue but absent from the single-cell reference (for example hair follicles or neutrophils)?")

CP(112, "Why do we use the `q05` cell abundance?")

READ("Save the annotated single-cell object", '''adata_sc.write_h5ad(f"{MY}/adata_sc_annotated_cell_state.h5ad")
''', ["Nothing to modify."])

CP(114, "From cell types to tissue: reading the Cell2location results")
READ("Load the Cell2location results into the spatial table", '''adata = sdata_all["table"]

q05 = pd.read_csv(Q05_FILE, index_col=0)          # spots x cell states, provided by the instructors
assert (q05.index == adata.obs_names).all(), "Spot names differ: use the spatial table of the instructors' checkpoint."

for state in q05.columns:
    adata.obs[f"c2l_{state}"] = q05[state].values          # one obs column per cell state
adata.obs["dominant_state"] = q05.idxmax(axis=1).astype("category")   # the most abundant state in each spot

adata = eh.relink_table(sdata_all, adata, SAMPLES)         # re-declare the table -> spots link
print(q05.shape)
q05.head()
''', ["Nothing to modify. Expected: a table with one row per spot and one column per cell state."])

READ("Save the spatial object with the Cell2location results", '''sdata_all.write(f"{MY}/spatial_data_annotated_with_cell2location.zarr", overwrite=True)
''', ["Nothing to modify (checkpoint, can take a minute)."])

# ============================================================================
# 6. READING THE RESULTS
# ============================================================================
CP(117, "We will draw many maps")
RUN("Helper: one variable on the four tissue sections", '''def spatial_row(color, axes, vmax=None, cmap="viridis"):
    """Draw one variable (an .obs column or a gene) on the 4 tissue sections.
    vmax: if given, all panels share the colour scale [0, vmax]."""
    eh.spatial_row(sdata_all, color, axes, SAMPLES, vmax=vmax, cmap=cmap)
''', ["Nothing to modify: from now on, spatial_row(variable, axes) draws a variable on the 4 sections."])

CP(120, "Step 1: does the deconvolution make sense?")
TURN("1a. How many cells per spot, and where?", """# What to do:
#   1. Store in adata.obs["c2l_total_cells"] the sum of q05 over the cell states of each spot.
#   2. Map it on the 4 sections with spatial_row (figure with 1 row x 4 columns).
#   3. Compute the Spearman correlation with the sequencing depth (adata.obs["total_counts"]) -> spearmanr.
# Expected: a map of cell density, and a printed correlation.

____
fig, axes = plt.subplots(1, len(SAMPLES), figsize=(5 * len(SAMPLES), 5))
____
fig.suptitle("Estimated number of cells per spot (sum of q05)")
fig.tight_layout()
plt.show()

____
print(f"Spearman(total estimated cells, total_counts) = {rho:.2f}")
""")
Q("Does the estimated cell density change in space and across samples? Is the correlation with sequencing depth high? Is that necessarily a problem?")

CP(122, "1b. Where is each cell state?")
TURN("1b. Where is each cell state? (shared colour scale)", """# What to do:
#   Choose 4-5 cell states that you find interesting (the list of available states is printed below),
#   and map each of them on the 4 sections, one row per state, using a colour scale shared by the
#   4 samples: vmax = 99th percentile of that state (q05[state].quantile(0.99)).
# Expected: a grid (states x samples). The plotting loop is written for you.

print(list(q05.columns))

states_to_show = [
    # "...",
]
assert len(states_to_show) > 0, "Choose at least one state"

fig, axes = plt.subplots(len(states_to_show), len(SAMPLES), figsize=(5 * len(SAMPLES), 5 * len(states_to_show)), squeeze=False)
for row, state in enumerate(states_to_show):
    vmax = q05[state].quantile(0.99)           # shared scale across the 4 samples
    spatial_row(f"c2l_{state}", axes[row], vmax=vmax)
    axes[row, 0].set_ylabel(state, fontsize=12)
fig.tight_layout()
plt.show()
""")
Q("Pick two states: where are they located, does it match the tissue (epidermis on top, dermis below, wound edge)? Do they appear or disappear over time? Did you see a state that looks uniformly spread everywhere, and what could that mean?")

MD("""
### 1c. An independent check with genes

This is the most convincing control. The signatures of the cell states come from **scRNA-seq**, but the gene expression we measured in each spot is an **independent measurement** made on the spatial data. If the deconvolution works, the abundance of a cell state must correlate with the expression of one of its marker genes **in the same spots**.

You will test a few (cell state, marker gene) pairs of your choice, and then look at the same comparison on the tissue: abundance on top, gene expression below. If the two maps overlap, the deconvolution is capturing real biology.
""")
TURN("1c. An independent check with genes", """# What to do:
#   Propose (cell state, marker gene) pairs, with genes NOT used to define the states (use the paper's table of
#   markers, or your own knowledge). For each pair, compute the Spearman correlation between the abundance
#   of the state (q05[state]) and the expression of the gene (adata.obs_vector(gene)) in the same spots.
#   One example is given. Add at least 3 more.
#   Then look at the same comparison on the tissue: abundance (top row) vs gene (bottom row).
# Expected: a table of correlations, and a figure with abundance and gene maps that overlap if the model works.

checks = [
    ("Keratinocytes_wound_activated", "KRT6A"),     # example given
    # ("...", "..."),
]

rows = []
for state, gene in checks:
    if gene not in adata.var_names:
        print(f"{gene} not in the filtered gene set, skipped")
        continue
    rho, p = spearmanr(q05[state].values, adata.obs_vector(gene))
    rows.append({"state": state, "marker_gene": gene, "spearman_rho": round(rho, 2)})
print(pd.DataFrame(rows).to_string(index=False))

# The same thing on the tissue, for the first pair: abundance (top) vs gene expression (bottom)
state, gene = checks[0]
fig, axes = plt.subplots(2, len(SAMPLES), figsize=(5 * len(SAMPLES), 10))
spatial_row(f"c2l_{state}", axes[0], vmax=q05[state].quantile(0.99))
spatial_row(gene, axes[1])
axes[0, 0].set_ylabel(f"Cell2location q05\\n({state})", fontsize=11)
axes[1, 0].set_ylabel(f"{gene} expression", fontsize=11)
fig.tight_layout()
plt.show()
""")
Q("Which pairs correlate well and which do not? For the ones that do not, can you think of a biological or technical reason?")

CP(126, "Step 2: linking the two annotations", replace=[
    ("our two **tentative labels**: cluster 3 (granulation tissue) and cluster 0 (inflamed tissue, Wound1). They now have", "the **tentative labels** you gave to the spatial clusters, especially the ones you were unsure about. They now have"),
    ("the name we gave them", "the name you gave them"),
])
RUN("Mean abundance of each cell state in each spatial cluster", """composition, concentration = eh.plot_cluster_by_state(q05, adata, cluster_key="leiden", name_key="spot_type")
""", ["Expected: two heatmaps (composition on the left, concentration on the right), one row per spatial cluster."])
Q("Does the cellular composition of each cluster confirm the name you gave it, or suggest another one? Focus on the clusters you were unsure about. (Remember: the reference contains no neutrophils.)")

CP(128, "Step 3: how does the composition change over time?")
RUN("Cell-state proportions per sample: spatial vs single-cell", """st_prop, sc_prop = eh.compare_state_proportions(q05, adata, adata_sc, SAMPLES)
""", ["Expected: the mean number of cells per spot per sample, then two stacked bar plots (ST left, scRNA-seq right)."])
Q("Describe how the composition changes from Skin to Wound30 in the two datasets. Where do they agree, where do they differ, and why?")

CP(130, "Step 4: who is next to whom?")
READ("Spatial graph for one sample", S(131), ["Nothing to modify. Change S to look at other samples later on."])
CP(132, "4a. Neighborhood enrichment")
TURN("4a. Neighborhood enrichment between spatial clusters", """# What to do:
#   Run the permutation test (sq.gr.nhood_enrichment, cluster_key="spot_type", seed=0) on the object `a`,
#   then plot the result (sq.pl.nhood_enrichment, same cluster_key, figsize=(7, 6)).
# Expected: a heatmap of z-scores: positive = neighbours more often than by chance, negative = avoidance.

____
____
""")
Q("Which pairs of spot types are strongly enriched as neighbours, and which avoid each other? Does it match the tissue architecture (epidermis / dermis)?")

CP(134, "4b. Spatial co-abundance of cell states")
READ("Spatial co-abundance of cell states", S(135), ["Nothing to modify. Look at how the hierarchical clustering groups the cell states."])
CP(136, "4c. How spatially structured")
READ("Moran's I of each cell state", S(137), ["Nothing to modify. Which states form precise structures and which are diffuse?"])
Q("Which cell states are the most and the least spatially structured? Is it what you expected from their biology?")

CP(138, "Take-home message and caveats",
   replace=[("- The ST cluster labels `spot_type` are working interpretations; clusters 0 and 3 are open questions that you are encouraged to explore.",
             "- Your `spot_type` labels are working interpretations: the deconvolution helps to test them, but it does not prove them.")])

# ----------------------------------------------------------------------------
# write the notebook
# ----------------------------------------------------------------------------
nb_s = {
    "cells": cells,
    "metadata": nb_t["metadata"],
    "nbformat": 4,
    "nbformat_minor": 5,
}
with open(STUDENT, "w") as f:
    json.dump(nb_s, f, indent=1, ensure_ascii=False)
    f.write("\n")

n_code = sum(c["cell_type"] == "code" for c in cells)
tags = {}
for c in cells:
    for t in c["metadata"]["tags"]:
        tags[t] = tags.get(t, 0) + 1
print(f"Written {STUDENT}: {len(cells)} cells ({n_code} code), tags: {tags}")
