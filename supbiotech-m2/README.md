# SupBioTech M2: scRNA-seq analysis of human skin wound healing

**Audience:** Master 2 students, SupBioTech (Paris)
**Data:** scRNA-seq only (GSE241132), one donor (PWH27) at D0 / D1 / D7 / D30
**Scope:** the whole scRNA-seq workflow, from the count matrices to cell types (first half), then a zoom on the keratinocytes: sub-states, pseudotime and gene sets (second half)

Biological background and data description: see the [main README](../README.md) and [`../data/README.md`](../data/README.md).

## Contents

| Path | What it contains | Published on the website |
|---|---|---|
| [`notebooks/Analysis_sc_students.ipynb`](notebooks/Analysis_sc_students.ipynb) | **Student notebook** (21 parts, two halves). Empty "Your turn" cells with a starter comment, and folded Hint 1 (the idea) / Hint 2 (full solution) | yes |
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
6. zoom in on a cell type (sub-clustering from raw counts), order cells along a pseudotime with an explicit root, and interpret gene lists with ORA and GSEA, knowing the traps of each.

## Outline

First half (the whole tissue): setup, 10x technology, `AnnData`, QC, doublets, normalisation, HVG, PCA, kNN graph and UMAP, batch effects, clustering (RESOLUTE), marker genes, annotation, validation, composition over time (Parts 1-16).
Second half (keratinocytes): zooming in, sub-clustering with and without Harmony, PAGA and diffusion pseudotime, gene trends, ORA and GSEA with GO Biological Process (Parts 17-21).

## Prerequisites and setup

```bash
mamba env create -f environment.yml
mamba activate supbiotech-sc
python -m ipykernel install --user --name skinwound --display-name "Python (skinwound)"
```

Data: see `../data/README.md` (or the course server path given by the instructor). Checkpoints of the instructor run go in `data/processed_data/SupBioTech/` (safety net for the students); the first run of Part 20 downloads and caches the GO library in `data/gene_sets/`.

## For the instructors

- The student notebook runs the same code as the instructor one (checked when it was built: replacing every "Your turn" cell by its Hint 2 gives exactly the instructor code). If you modify one of the two notebooks, update the other by hand.
- The instructor run must exist before the course: `adata_sc_leiden_res0.2.h5ad`, `adata_sc_leiden_res0.2_annotated.h5ad`, `adata_sc_annotated_cell_state.h5ad` and `sc_data_PWH27_processed.h5ad` in `processed_data/SupBioTech/`.
- Cluster numbers are not reproducible across machines: the reference solutions translate them through the marker genes.

## Instructors
<!-- TODO -->
