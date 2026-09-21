# Data

Raw data are **not** versioned (see `.gitignore`); they are public on GEO.

## Downloads
- **Visium (spatial)**: GEO series [GSE241124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241124) → *Supplementary files* (`GSE241124_RAW`, plus `GSE241124_spatialseq_metadata_acutewound.txt`).
- **Chromium (scRNA-seq)**: GEO series [GSE241132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241132) → *Supplementary files* (one `.zip` per sample, plus `GSE241132_cell_metadata.txt`).

Expected layout after download:

```
data/
├── GSE241124_RAW/
│   ├── GSE241124_spatialseq_metadata_acutewound.txt
│   └── <sample>/filtered_feature_bc_matrix.h5  +  <sample>/spatial/…
└── GSE241132_RAW/
    ├── GSE241132_cell_metadata.txt
    └── GSM77170xx_PWHxxDx.zip   (12 samples)
```

> On the shared course server the data are already available at: `<TODO: path>`.

## Visium samples (from `GSE241124_spatialseq_metadata_acutewound.txt`)

| Sample id (`orig.ident`) | Donor | Patient | Condition | Sex | Age | Seq. batch | Spots (after authors' QC) |
|---|---|---|---|---|---|---|---|
| P17401_1001 | Donor1 | PWH19 | Skin | M | 24 | 1 | 810 |
| P17401_1002 | Donor1 | PWH19 | Wound1 | M | 24 | 1 | 1832 |
| P17401_1003 | Donor1 | PWH19 | Wound7 | M | 24 | 1 | 1689 |
| P17401_1004 | Donor1 | PWH19 | Wound30 | M | 24 | 1 | 1278 |
| P20063_101 | Donor2 | PWH20 | Skin | F | 46 | 2 | 842 |
| P20063_102 | Donor2 | PWH20 | Wound1 | F | 46 | 2 | 1758 |
| P20063_103 | Donor2 | PWH20 | Wound7 | F | 46 | 2 | 1606 |
| P20063_104 | Donor2 | PWH20 | Wound30 | F | 46 | 2 | 2205 |
| P20063_105 | Donor3 | PWH28 | Skin | M | 24 | 3 | 685 |
| P20063_106 | Donor3 | PWH28 | Wound1 | M | 24 | 3 | 1317 |
| P20063_107 | Donor3 | PWH28 | Wound7 | M | 24 | 3 | 1537 |
| P20063_108 | Donor3 | PWH28 | Wound30 | M | 24 | 3 | 2019 |
| P20063_109 | Donor4 | PWH27 | Skin | F | 22 | 3 | 740 |
| P20063_110 | Donor4 | PWH27 | Wound1 | F | 22 | 3 | 1608 |
| P20063_111 | Donor4 | PWH27 | Wound7 | F | 22 | 3 | 1212 |
| P20063_112 | Donor4 | PWH27 | Wound30 | F | 22 | 3 | 1777 |

Total: 22,915 spots (matches the paper). The metadata file also holds the authors' spot annotation (`AnnoType`).

> **TODO (to verify):** in the downloaded folders, samples `P20063_109`–`P20063_112` appear under the names `Skin`, `Wound1`, `Wound7`, `Wound30`. Confirm this correspondence (e.g. by comparing spot barcodes with the metadata) before teaching it.

## scRNA-seq samples (from `GSE241132_cell_metadata.txt`)

| Patient | Skin (D0) | Wound1 (D1) | Wound7 (D7) | Wound30 (D30) |
|---|---|---|---|---|
| PWH26 | PWH26D0 (4,242 cells) | PWH26D1 (4,709) | PWH26D7 (4,483) | PWH26D30 (3,408) |
| PWH27 | PWH27D0 (3,265) | PWH27D1 (6,491) | PWH27D7 (6,518) | PWH27D30 (5,953) |
| PWH28 | PWH28D0 (4,083) | PWH28D1 (4,578) | PWH28D7 (5,358) | PWH28D30 (5,735) |

Total: 58,823 cells after the authors' QC (matches the paper). Cell metadata include the authors' annotation (`newCellTypes`, `newMainCellTypes`), useful as a reference to compare with the students' own annotation.

Note that **PWH28 (Donor3) and PWH27 (Donor4) appear in both technologies**, which is what allows the ST + scRNA-seq integration.
