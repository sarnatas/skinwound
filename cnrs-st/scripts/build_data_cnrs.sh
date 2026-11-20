#!/usr/bin/env bash
# Builds data/DATA_CNRS (clean names) by COPYING from the raw GEO downloads. Originals are never touched.
# Usage: bash cnrs-st/scripts/build_data_cnrs.sh   (from the repo root)
set -euo pipefail
ROOT="$(pwd)"; SRC_ST="$ROOT/data/GSE241124_RAW"; SRC_SC="$ROOT/data/GSE241132_RAW"; DST="$ROOT/data/DATA_CNRS"
mkdir -p "$DST/visium" "$DST/scrna"

# --- 1. Visium: GEO folder -> DonorN_condition ---
declare -A MAP=( [P17401_1001]=Donor1_Skin [P17401_1002]=Donor1_Wound1 [P17401_1003]=Donor1_Wound7 [P17401_1004]=Donor1_Wound30
                 [P20063_101]=Donor2_Skin  [P20063_102]=Donor2_Wound1  [P20063_103]=Donor2_Wound7  [P20063_104]=Donor2_Wound30
                 [P20063_105]=Donor3_Skin  [P20063_106]=Donor3_Wound1  [P20063_107]=Donor3_Wound7  [P20063_108]=Donor3_Wound30
                 [Skin]=Donor4_Skin        [Wound1]=Donor4_Wound1      [Wound7]=Donor4_Wound7      [Wound30]=Donor4_Wound30 )
for k in "${!MAP[@]}"; do
  d="$DST/visium/${MAP[$k]}"; mkdir -p "$d/spatial"
  cp "$SRC_ST/$k/filtered_feature_bc_matrix.h5" "$d/"
  sp=$(find "$SRC_ST/$k" -mindepth 1 -maxdepth 1 -type d ! -name __MACOSX | head -1)
  cp "$sp"/{tissue_hires_image.png,tissue_lowres_image.png,tissue_positions_list.csv,scalefactors_json.json} "$d/spatial/"
done
cp "$SRC_ST/GSE241124_spatialseq_metadata_acutewound.txt" "$DST/visium/"

# --- 2. scRNA-seq (GSE241132): PWH27 + PWH28 (the donors also present in Visium), 10x names matrix/barcodes/features ---
tmp=$(mktemp -d)
for p in PWH27 PWH28; do for t in D0 D1 D7 D30; do
  z=$(ls "$SRC_SC"/GSM*_${p}${t}.zip); unzip -q -o "$z" -x '__MACOSX/*' -d "$tmp"
  o="$DST/scrna/${p}${t}"; mkdir -p "$o"
  for f in matrix.mtx barcodes.tsv features.tsv; do mv "$tmp/${p}${t}/${p}${t}_${f}.gz" "$o/${f}.gz"; done
done; done
rm -rf "$tmp"
cp "$SRC_SC/GSE241132_cell_metadata.txt" "$DST/scrna/"

# --- 3. Paper marker genes ---
cp "$ROOT/data/paper_marker_genes.csv" "$DST/"
echo "Done: $DST"; du -sh "$DST"
