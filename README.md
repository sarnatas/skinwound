# Spatiotemporal roadmap of human skin wound healing
### Teaching material built on published Visium (spatial) and Chromium (scRNA-seq) datasets

This repository contains practical courses on **spatial transcriptomics and single-cell RNA-seq**, built on the dataset of

> Liu Z., Bian X., Luo L., ..., Sommar P., Li D., Xu Landén N. **Spatiotemporal single-cell roadmap of human skin wound healing.** *Cell Stem Cell* 32, 479–498 (2025). doi: [10.1016/j.stem.2024.11.013](https://doi.org/10.1016/j.stem.2024.11.013)

The authors created acute wounds in healthy volunteers and sampled the wound edge before wounding (`Skin`) and at day 1, 7 and 30 (`Wound1`, `Wound7`, `Wound30`), profiling the same individuals with 10x Visium and 10x Chromium scRNA-seq. The time course covers the inflammatory, proliferative and remodeling phases of human wound healing. The courses use these public data to teach quality control, clustering, batch effects, cell-type annotation, spatial deconvolution (Cell2location) and spatial neighbourhood analysis, with an emphasis on critical thinking about analysis choices.

## Practical notebook (ENS de Lyon, NGS practicals)

**[Open the notebook as a website](https://sarnatas.github.io/skinwound/ens-ngs/notebooks/Analysis_open.html)**: navigable, with table of contents (instructions and hints; results appear when you run it).

- Course page, setup and data layout: <https://sarnatas.github.io/skinwound/ens-ngs/>
- To work on it: download [`Analysis_open.ipynb`](ens-ngs/notebooks/Analysis_open.ipynb) (a single file) and create the environment (see below).

The notebook is a guide: students explore the tools themselves, with folded hints for the technically difficult steps, and deliver their own notebook at the end. The website is rebuilt automatically from it at every push to `main`.

## Courses

| Course | Audience | Data used | Folder |
|---|---|---|---|
| NGS practicals | ENS de Lyon | Spatial transcriptomics **and** its integration with scRNA-seq | [`ens-ngs/`](ens-ngs/) |
| Master 2 course | SupBioTech, Paris | scRNA-seq | [`supbiotech-m2/`](supbiotech-m2/) *(in preparation)* |
| Professional training | CNRS formation entreprise | Spatial transcriptomics only | [`cnrs-st/`](cnrs-st/) *(in preparation)* |

## Data

| Accession | Technology | Content | Used by |
|---|---|---|---|
| [GSE241124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124) | 10x Visium | Human acute wounds, 4 donors × 4 conditions (16 sections, 22,915 spots) | `ens-ngs`, `cnrs-st` |
| [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132) | 10x Chromium 3' v3.1 | Human acute wounds, 3 donors × 4 conditions (12 samples, 58,823 cells) | `ens-ngs`, `supbiotech-m2` |

Raw data are **not** stored in this repository; download instructions and sample tables are in [`data/README.md`](data/README.md). Biology of wound healing, experimental design and the main findings of the paper are summarised on the [Background page](https://sarnatas.github.io/skinwound/background.html).

## Repository layout

```
.
├── README.md            <- you are here
├── _quarto.yml          <- configuration of the website (Quarto)
├── index.qmd, background.qmd, styles.css   <- pages of the website
├── environment.yml      <- shared Python environment
├── data/                <- download instructions only (data are git-ignored)
├── common/              <- code shared by the courses
├── figures/             <- figures reused in several courses
├── ens-ngs/             <- ST + scRNA-seq integration   (notebooks/, slides/, scripts/, index.qmd)
├── supbiotech-m2/       <- scRNA-seq
└── cnrs-st/             <- ST only
```

## Software environment

```bash
mamba env create -f environment.yml
mamba activate skinwound
```

Main packages: scanpy, squidpy, spatialdata (+ spatialdata-io, spatialdata-plot), anndata, leidenalg, Cell2location. See `environment.yml` for exact versions.

## The website

The site is built with [Quarto](https://quarto.org) and published on GitHub Pages by the workflow `.github/workflows/publish.yml`. Only the pages listed in `_quarto.yml` are published (the practical notebook, never the instructor one). Notebooks are rendered, **not executed**. To preview it locally:

```bash
quarto preview
```

## How to cite / reuse

Please cite the original paper and the GEO accessions when using these data. Course material: see `LICENSE` <!-- TODO: choose a license, e.g. CC BY 4.0 for the material and MIT for the code -->.

**Author of the courses:** Sergio Sarnataro <!-- TODO: affiliation / contact -->
