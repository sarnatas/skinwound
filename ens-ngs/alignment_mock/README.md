# Mock alignment with Space Ranger (toy data)

A short warm-up: run `spaceranger count` on **tiny synthetic FASTQ files** for the four
samples of the practical (Skin, Wound1, Wound7, Wound30). The goal is to *see* what an
alignment run looks like (inputs, command line, logs, QC report), **not** to produce
results. Afterwards, forget the output and go back to the notebook, which uses the real,
already aligned data (GSE241124).

```bash
cd ens-ngs/alignment_mock
bash alignment_mock.sh
```

The script is just a loop over the four samples calling `spaceranger count`; read it
first, it is ~15 lines. The toy reference `toy_data/toy_ref` is assumed to exist already (built
once by the instructor, see the commented `mkref` command at the top of the script).
Adapt `--localcores` / `--localmem` to your machine.

## What is in this folder

| Path | What it is |
|---|---|
| `alignment_mock.sh` | the script students run |
| `toy_data/fastq/<sample>/` | `<sample>_S1_L001_R1_001.fastq.gz` / `R2` (~42k read pairs, ~1.5 MB each) |
| `toy_data/images/<sample>.jpg` | tissue images (JPEG, 2000 px; Space Ranger does not take PNG) |
| `toy_data/reference_src/` | FASTA + GTF of the toy genome: only the *input* of `mkref`, not a valid `--transcriptome` |
| `toy_data/toy_ref/` | the reference built by `spaceranger mkref` (instructor, once); this is what `--transcriptome` points to |
| `make_toy_data.py` | instructor-side generator (deterministic). Students do not need it |

## About the toy data

* **Real:** spot barcodes (from each sample's `tissue_positions_list.csv`) and the
  relative expression of 47 skin/wound marker genes per spot (thinned to ~25k UMIs per
  sample).
* **Synthetic:** every nucleotide. Gene names are real (KRT5, S100A8, COL1A1…) but the
  sequences are random: this is **not** the human genome. Do not read anything biological
  into the toy results.
* Visium v1 layout: R1 = 16 nt barcode + 12 nt UMI, R2 = 90 nt cDNA, with PCR duplicates,
  barcode errors, ambient RNA and unmappable reads, so the QC report has something to show.

## What to look at

Open `<sample>/outs/web_summary.html`: mapping rate, valid barcodes, sequencing
saturation, spots under tissue, median genes per spot (≤ 47 by construction).

## Status

The reads and reference were validated with STARsolo (97% uniquely mapped reads; recovered
UMI matrix matches the simulated truth). The script itself has **not yet been run with the
real Space Ranger**. Things to check on the first run: automatic fiducial detection on the
2000-px GEO image (converted to JPEG) (`--unknown-slide=visium-1`; fallback: `--loupe-alignment=file.json`)
and run time.
