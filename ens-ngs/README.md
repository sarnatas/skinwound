# ENS – NGS practicals: spatial transcriptomics and integration with scRNA-seq

**Audience:** Students of the NGS practicals, ENS de Lyon
**Data:** Visium (GSE241124) **and** scRNA-seq (GSE241132)
**Scope:** Analysis of the Visium data, then integration with the scRNA-seq analysis (cell-type mapping / deconvolution)

Biological background and data description: see the [main README](../README.md) and [`../data/README.md`](../data/README.md).

## Contents

| Folder | What it contains | Status |
|---|---|---|
| [`slides/`](slides/) | Lecture slides (PDF / HTML) | TODO |
| [`notebooks/`](notebooks/) | Analysis notebooks (marimo `.py` and/or Jupyter `.ipynb`) | TODO |
| [`site/`](site/) | Static website (rendered notebooks, HTML exports) served by GitHub Pages | placeholder |

Online version: `https://<github-user>.github.io/<repo-name>/ens-ngs/site/`  <!-- TODO -->

## Learning objectives
<!-- TODO: 4-6 bullet points -->

## Planned outline
1. Loading Visium data (spatialdata-io / scanpy) and quality control
2. Normalisation, dimensionality reduction, clustering of spots
3. Spatial visualisation and spatially variable genes (squidpy)
4. Reference scRNA-seq: annotation of cell types
5. Integration ST + scRNA-seq: mapping / deconvolution of cell types onto spots
6. Spatial neighbourhoods and interpretation across healing phases (Skin, Wound1, Wound7, Wound30)

## Prerequisites and setup
Python environment: see the root `environment.yml`.
```bash
mamba env create -f ../environment.yml
mamba activate skinwound
```
Data: see `../data/README.md` (or the course server path given by the instructor).

## Instructors
<!-- TODO -->
