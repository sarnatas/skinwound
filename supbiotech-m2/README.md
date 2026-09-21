# SupBioTech M2: scRNA-seq analysis of human skin wound healing

**Audience:** Master 2 students, SupBioTech (Paris)
**Data:** scRNA-seq only (GSE241132)
**Scope:** Full scRNA-seq workflow, from the count matrices to cell types and differences between healing phases

Biological background and data description: see the [main README](../README.md) and [`../data/README.md`](../data/README.md).

## Contents

| Folder | What it contains | Status |
|---|---|---|
| [`slides/`](slides/) | Lecture slides (PDF / HTML) | TODO |
| [`notebooks/`](notebooks/) | Analysis notebooks (marimo `.py` and/or Jupyter `.ipynb`) | TODO |
| [`site/`](site/) | Static website (rendered notebooks, HTML exports) served by GitHub Pages | placeholder |

Online version: `https://<github-user>.github.io/<repo-name>/supbiotech-m2/site/`  <!-- TODO -->

## Learning objectives
<!-- TODO: 4-6 bullet points -->

## Planned outline
1. Loading 10x matrices and quality control (genes/cell, UMIs, % mitochondrial)
2. Filtering, normalisation, highly variable genes
3. PCA, neighbours, UMAP, Leiden clustering
4. Batch effects across donors and integration
5. Marker genes and cell-type annotation (compare with the authors' annotation)
6. Differential abundance and expression across Skin / Wound1 / Wound7 / Wound30

## Prerequisites and setup
Python environment: see the root `environment.yml`.
```bash
mamba env create -f ../environment.yml
mamba activate skinwound
```
Data: see `../data/README.md` (or the course server path given by the instructor).

## Instructors
<!-- TODO -->
