# Spatiotemporal roadmap of human skin wound healing
### Teaching material built on the published Visium (ST) and Chromium (scRNA-seq) datasets

This repository gathers the material of three courses that use the same published dataset of human acute skin wound healing:

| Course | Audience | Data used | Folder |
|---|---|---|---|
| NGS practicals | ENS de Lyon | Spatial transcriptomics **and** its integration with scRNA-seq | [`ens-ngs/`](ens-ngs/) |
| Master 2 course | SupBioTech, Paris | scRNA-seq | [`supbiotech-m2/`](supbiotech-m2/) |
| Professional training | CNRS formation entreprise | Spatial transcriptomics only | [`cnrs-st/`](cnrs-st/) |

Each course folder is self-contained (slides, static website, notebooks). This page only covers what all three share: the biology and the data.
If you follow one course, start from the README of your folder.

> Site (GitHub Pages): `https://<github-user>.github.io/<repo-name>/`  <!-- TODO: fill in -->

---

## 1. The source study

Liu Z., Bian X., Luo L., ..., Sommar P., Li D., Xu Landén N. **Spatiotemporal single-cell roadmap of human skin wound healing.**
*Cell Stem Cell* 32, 479–498 (2025). doi: [10.1016/j.stem.2024.11.013](https://doi.org/10.1016/j.stem.2024.11.013) (open access, CC BY).

The authors combined scRNA-seq (10x Chromium) and spatial transcriptomics (10x Visium) on human skin wounds created *in vivo* in healthy volunteers, and compared them with chronic wounds (venous and diabetic foot ulcers) and mouse wounds.
Resources released by the authors: interactive atlas at <https://www.xulandenlab.com/tools>, original analysis code at Zenodo [10.5281/zenodo.14176654](https://doi.org/10.5281/zenodo.14176654).

## 2. Biological background

### Wound healing in three overlapping phases
1. **Inflammation** (first days): recruitment of neutrophils and monocytes/macrophages, clearing of debris and microbes.
2. **Proliferation** (about one week): re-epithelialization (keratinocytes migrate and proliferate to cover the wound), formation of granulation tissue (fibroblasts, new capillaries).
3. **Remodeling** (weeks to months): extracellular matrix reorganization, scar maturation.

Failed re-epithelialization is the hallmark of chronic wounds such as diabetic foot ulcers (DFU) and venous ulcers (VU). Rodent models are widely used but differ markedly from human skin (thicker human epidermis, fewer hair follicles, no *panniculus carnosus*-driven contraction in humans), hence the value of a human time course.

### Experimental design
- Healthy volunteers (Karolinska University Hospital, Stockholm). Three full-thickness wounds were made on the upper buttock of each donor with a 4 mm punch biopsy.
- The wound edge was sampled with a 6 mm punch at **day 1, day 7 and day 30** after injury, plus **intact skin** before wounding, all from the *same individuals*. The four conditions are named `Skin`, `Wound1`, `Wound7`, `Wound30` and stand for the inflammatory, proliferative and remodeling phases (the paper supports this assignment with GO analysis of DEGs and comparison with bulk RNA-seq).
- scRNA-seq: 3 donors. Visium ST: 4 donors. Two donors were analysed with both technologies.
- Low donor heterogeneity was observed: samples cluster by condition rather than by donor or by technology (pseudobulk PCA, without batch correction), which makes the dataset well suited to compare healing stages.

### Main findings (useful hooks for the exercises)
- **Cell atlas.** 58,823 cells (12 samples) in 27 clusters and 9 main cell types: keratinocytes, fibroblasts, myeloid and lymphoid cells, endothelial cells, mast cells, pericytes/smooth muscle, melanocytes, Schwann cells. Canonical markers include KRT5/KRT10 (keratinocytes), COL1A1 (fibroblasts), PECAM1 (endothelium), LYZ/HLA-DRA (myeloid), CD3D/NKG7 (lymphoid).
- **Spatial map.** 22,915 Visium spots (16 sections), 17 spatial clusters: basal and suprabasal epidermis, hair follicle, a **wound-edge cluster** (KRT6A/B/C, KRT16, KRT17, S100A8/9), fibroblast-rich dermis, sweat and sebaceous glands, immune, endothelial, smooth muscle and mast-cell clusters. Cell types from the scRNA-seq were mapped onto the spots by deconvolution (cell2location), then grouped into 15 spatial "niches" with NMF.
- **Re-epithelialization.** A non-proliferative migrating front (basal-migrating Bas-mig and spinous-migrating Spi-mig keratinocytes, expressing MMP1/MMP3, KRT6, S100A proteins) is surrounded by a proliferating hub (Bas-prolif). Unlike mouse, human wounds show a clear separation between migration and proliferation.
- **FOSL1** (AP-1 component) is identified as a master regulator of keratinocyte migration (SCENIC regulons, in silico perturbation with CellOracle, siRNA scratch assay).
- **A relay race.** Pro-inflammatory macrophages (Mac_inf, Mac1; CXCL1/CXCL5 and EREG signals) support keratinocyte migration at the inflammatory stage, then proliferating fibroblasts (HGF, FGF2, TGFB1) take over at the proliferative stage.
- **Chronic wounds.** Migrating keratinocytes are strongly reduced in DFU and absent in VU, and inflammatory responses are impaired rather than simply excessive.
- **Human vs mouse.** Conserved migration markers (FOSL1, AREG, IL24, NRG1, GJB2) but many human-specific genes (MMP1, S100A2/7/8/9, SERPINB3/4).

### Methodological limitation to keep in mind
Visium spots (55 µm) are **not single cells**: each spot mixes several cells (the authors assumed about 20 cells per spot for deconvolution). This is exactly why integrating ST with scRNA-seq is instructive.

## 3. Datasets

| Accession | Technology | Content | Used by |
|---|---|---|---|
| [GSE241124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124) | 10x Visium | Human acute wounds, 4 donors × 4 conditions (16 sections) | `ens-ngs`, `cnrs-st` |
| [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132) | 10x Chromium 3' v3.1 | Human acute wounds, 3 donors × 4 conditions (12 samples) | `ens-ngs`, `supbiotech-m2` |

Other datasets of the study (not used in the courses): GSE265972 (venous ulcer scRNA-seq), GSE218430 (mouse wounds).
Downloading and describing the files is explained in [`data/README.md`](data/README.md). Raw data are **not** stored in this repository.

## 4. Repository layout

```
.
├── README.md            <- you are here (biology + data)
├── index.html           <- landing page of the GitHub Pages site
├── environment.yml      <- shared Python environment
├── data/                <- download instructions only (data are git-ignored)
├── common/              <- code shared by the courses
├── figures/             <- figures reused in several courses
├── ens-ngs/             <- ST + scRNA-seq integration   (notebooks/, slides/, site/)
├── supbiotech-m2/       <- scRNA-seq                    (notebooks/, slides/, site/)
└── cnrs-st/             <- ST only                      (notebooks/, slides/, site/)
```

## 5. Software environment

```bash
mamba env create -f environment.yml
mamba activate skinwound
```

Main packages: scanpy, squidpy, spatialdata (+ spatialdata-io, spatialdata-plot), anndata, leidenalg, and [marimo](https://marimo.io) as reactive notebook. See `environment.yml` for exact versions.

## 6. How to cite / reuse

Please cite the original paper (above) and the GEO accessions when using these data. Course material: see `LICENSE` <!-- TODO: choose a license, e.g. CC BY 4.0 for the material and MIT for the code -->.

**Author of the courses:** Sergio Sarnataro <!-- TODO: affiliation / contact -->
