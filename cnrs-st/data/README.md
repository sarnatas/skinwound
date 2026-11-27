# Data for the CNRS course: `sc_reference_mod.h5ad`

`sc_reference_mod.h5ad` is the scRNA-seq reference used in Part D4 (sc–st integration: signature scores and Tangram) of the
reference notebook [`../notebooks/CNRS_ST_reference.ipynb`](../notebooks/CNRS_ST_reference.ipynb).

**Download:** [sc_reference_mod.h5ad](https://github.com/sarnatas/skinwound/raw/main/cnrs-st/data/sc_reference_mod.h5ad) (66 MB).
Put it in your `DATA` folder; the notebook reads it with `sc.read_h5ad(f"{DATA}/sc_reference_mod.h5ad")`.

## Source and references

This file is a **subset of public data**, redistributed here for teaching only. If you use these data, cite the original study:

> Liu Z, Bian X, Luo L, et al. *Spatiotemporal single-cell roadmap of human skin wound healing.*
> **Cell Stem Cell** 32, 479–498 (2025). <https://doi.org/10.1016/j.stem.2024.11.013>
> (open access article under the CC BY 4.0 license)

- Raw count matrices and cell metadata (authors' annotation): GEO [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132) (10x Chromium scRNA-seq).
- The Visium data analysed in the notebook are not stored here: GEO [GSE241124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124).

## Content

`AnnData`, 7,710 cells × 36,600 genes.

| Slot | Content |
|---|---|
| `X` | raw UMI counts (not normalised) |
| `obs["Patient"]` | donor: PWH27 (= Visium Donor4) or PWH28 (= Visium Donor3), the two donors present in both technologies |
| `obs["Condition"]` | time point (healthy skin, wound day 1, day 7, day 30) |
| `obs["newCellTypes"]` | the authors' fine annotation (cell states) |
| `obs["cell_type"]` | the authors' main cell types (`newMainCellTypes` in GEO), with names without spaces or `&` (e.g. `Pericyte_Smooth_muscle`, `Hair_follicle`) |

Cell counts per `cell_type`: Keratinocyte 2,400 · Myeloid 1,317 · Fibroblast 1,200 · Lymphoid 900 · Endothelial 547 · Hair_follicle 300 · Mast 300 · Melanocyte 300 · Pericyte_Smooth_muscle 300 · Schwann 146.

## How it was built

With [`../scripts/build_sc_reference.py`](../scripts/build_sc_reference.py) from the GEO downloads:

1. cells of PWH27 and PWH28, all four time points (`GSE241132_cell_metadata.txt`, authors' QC);
2. at most 300 cells per fine label (`newCellTypes`), random sampling with `random_state=0`;
3. raw counts read from the 10x matrices of the 8 samples, cell barcodes prefixed with the sample name (`PWH27D0_AAAC...`);
4. annotation joined from the authors' metadata; `newMainCellTypes` then renamed `cell_type`, with spaces and `&` replaced by `_`
   (names with spaces or `&` are refused by `SpatialData`).
