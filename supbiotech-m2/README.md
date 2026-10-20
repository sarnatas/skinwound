# SupBioTech M2: scRNA-seq analysis of human skin wound healing

**Audience:** Master 2 students, SupBioTech (Paris)
**Data:** scRNA-seq only (GSE241132), one donor (PWH27) at D0 / D1 / D7 / D30
**Scope:** the whole scRNA-seq workflow, from the count matrices to cell types (first half), then a zoom on the keratinocytes: sub-states and pseudotime (second half)

Biological background and data description: see the [main README](../README.md) and [`../data/README.md`](../data/README.md).

## Contents

| Path | What it contains | Published on the website |
|---|---|---|
| [`notebooks/Analysis_sc_students.ipynb`](notebooks/Analysis_sc_students.ipynb) | **Student notebook** (20 parts, two halves). Empty "Your turn" cells with a starter comment, and folded Hint 1 (the idea) / Hint 2 (full solution) | yes |
| `notebooks/Analysis_sc_full.ipynb` | **Instructor notebook**: same content, all code visible, no outputs. Single runnable file = basic + advanced analysis | no |
| `notebooks/Analysis_sc.ipynb`, `Analysis_sc_advanced.ipynb`, `Analysis_sc_advanced_students.ipynb` | the earlier separate notebooks (first half with outputs, second half teacher / student). Superseded by the two files above | no |
| [`environment.yml`](environment.yml) | Python environment of the course (`supbiotech-sc`) | yes |
| [`slides/`](slides/) | Lecture slides | no |
| [`index.qmd`](index.qmd) | Course page of the website (setup, data layout, how we work) | yes |

Online version: https://sarnatas.github.io/skinwound/supbiotech-m2/index.html

## Learning objectives

By the end of the course the student can:

1. explain what a scRNA-seq count matrix is and what each slot of an `AnnData` object contains;
2. choose and justify QC thresholds, and remove doublets and uninformative genes;
3. explain what normalisation, highly variable genes, PCA, the kNN graph and UMAP do, and what a UMAP can *not* tell;
4. decide whether a batch effect should be corrected, given what the "batch" variable represents;
5. cluster the cells, annotate the clusters from marker genes, and validate the annotation against an independent reference;
6. zoom in on a cell type (sub-clustering from raw counts), and order cells along a pseudotime with an explicit root.

## Outline

First half (the whole tissue): setup, 10x technology, `AnnData`, QC, doublets, normalisation, HVG, PCA, kNN graph and UMAP, batch effects, clustering (RESOLUTE), marker genes, annotation, validation, composition over time (Parts 1-16).
Second half (keratinocytes): zooming in, sub-clustering with and without Harmony, PAGA and diffusion pseudotime, gene trends (Parts 17-20).

## Prerequisites and setup

```bash
conda env create -f environment.yml
conda activate supbiotech-sc
python -m ipykernel install --user --name skinwound --display-name "Python (skinwound)"
```

Data: the students download `DATA.zip` (course archive, ~300 MB: the four samples `PWH27D0/D1/D7/D30`, `GSE241132_cell_metadata.txt` and `paper_marker_genes.csv`, directly inside the `DATA/` folder, no `GSE241132_RAW` sub-folder) and unzip it. The notebook has two paths to set in the setup cell (`DATA` for the input data, `OUT` for the results) and checks them before going on. Students install conda with **Miniforge** (open source) and work locally; step-by-step instructions for Windows, Mac and Linux are on the course page ([`index.qmd`](index.qmd)).

Platforms: the pinned `environment.yml` was tested on Linux; for Windows and Apple-silicon Macs the availability of every package was checked, not a full run. On **Windows** `harmonypy 0.0.10` is used (version 2.x has no Windows package). **Macs with an Intel processor cannot install it** (no `numba` package for it): those students work in a group.

## For the instructors

- The student notebook runs the same code as the instructor one (checked when it was built: replacing every "Your turn" cell by its Hint 2 gives exactly the instructor code). If you modify one of the two notebooks, update the other by hand.
- Cluster numbers are not reproducible across machines: the reference solutions translate them through the marker genes.

## Instructors
<!-- TODO -->
