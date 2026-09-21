# CNRS professional training: spatial transcriptomics (Visium)

**Audience:** Professionals (CNRS formation entreprise)
**Data:** Visium only (GSE241124)
**Scope:** Spatial transcriptomics only: no scRNA-seq required

Biological background and data description: see the [main README](../README.md) and [`../data/README.md`](../data/README.md).

## Contents

| Folder | What it contains | Status |
|---|---|---|
| [`slides/`](slides/) | Lecture slides (PDF / HTML) | TODO |
| [`notebooks/`](notebooks/) | Analysis notebooks (marimo `.py` and/or Jupyter `.ipynb`) | TODO |
| [`site/`](site/) | Static website (rendered notebooks, HTML exports) served by GitHub Pages | placeholder |

Online version: `https://<github-user>.github.io/<repo-name>/cnrs-st/site/`  <!-- TODO -->

## Learning objectives
<!-- TODO: 4-6 bullet points -->

## Planned outline
1. Principles of Visium and structure of the data (matrix + image + spot coordinates)
2. Loading, quality control and filtering of spots
3. Normalisation, clustering, spatial visualisation
4. Spatially variable genes and marker genes of tissue regions (epidermis, dermis, wound edge)
5. Comparison across healing phases
6. Interpretation and limits of spot-based resolution

## Prerequisites and setup
Python environment: see the root `environment.yml`.
```bash
mamba env create -f ../environment.yml
mamba activate skinwound
```
Data: see `../data/README.md` (or the course server path given by the instructor).

## Instructors
<!-- TODO -->
