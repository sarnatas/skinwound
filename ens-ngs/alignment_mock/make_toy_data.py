#!/usr/bin/env python3
"""
make_toy_data.py -- generate the TOY reference + TOY FASTQ files used by
alignment_mock.sh (instructor-side script: students do NOT need to run it).

What is generated (all deterministic: same seed -> byte-identical files)
-----------------------------------------------------------------------
toy_data/reference_src/toy_genome.fa synthetic genome (2 contigs, ~45 "genes")
toy_data/reference_src/toy_genes.gtf matching annotation
toy_data/fastq/<sample>/<sample>_S1_L001_R{1,2}_001.fastq.gz
                                   Visium v1 style read pairs
                                   R1 = 16 nt spot barcode + 12 nt UMI (28 nt)
                                   R2 = 90 nt cDNA (sense strand, 3' biased)

What is REAL and what is SYNTHETIC
----------------------------------
* REAL     : the spot barcodes (taken from each sample's tissue_positions_list.csv,
             i.e. the Visium v1 whitelist) and the *relative* expression of the
             panel genes per spot (taken from filtered_feature_bc_matrix.h5 and
             randomly thinned to a few thousand UMIs, so the toy data keep a
             recognisable spatial pattern).
* SYNTHETIC: every nucleotide of the genome/reads. Gene *names* are real
             (skin / wound-healing markers) but their sequences are random, so
             the toy reference is NOT the human genome.

Usage
-----
    python make_toy_data.py --data-dir ../../data/GSE241124_RAW
    python make_toy_data.py --data-dir ... --target-umis 40000   # bigger files

Requires only numpy, pandas and h5py (all in the course environment).
"""
import argparse
import gzip
import io
import os
import random
import sys

import h5py
import numpy as np
import pandas as pd

SAMPLES = ["Skin", "Wound1", "Wound7", "Wound30"]

# Panel of genes present in the toy reference. Names come from
# data/paper_marker_genes.csv (major cell types, ST categories) plus the
# wound-edge markers and three housekeeping genes. Sequences are random.
PANEL = [
    "KRT5", "KRT14", "KRT15", "COL17A1", "KRT1", "KRT10", "KRT2", "LORICRIN",
    "KRT6A", "KRT6B", "KRT16", "KRT17", "S100A8", "S100A9", "KRT75", "KRT77",
    "DCD", "FADS2", "KRT79", "COL1A1", "COL1A2", "POSTN", "MMP2", "FBLN1",
    "ADAM12", "LYZ", "HLA-DRA", "CD3D", "NKG7", "FABP4", "CD163", "PECAM1",
    "VWF", "LYVE1", "TPSAB1", "TPSB2", "ACTA2", "TAGLN", "MYH11", "MYL9",
    "TYRP1", "PMEL", "SOX10", "MKI67", "ACTB", "GAPDH", "B2M",
]

READ1_LEN = 28      # 16 bp barcode + 12 bp UMI
BC_LEN = 16
UMI_LEN = 12
READ2_LEN = 90
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n, gc=0.45):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]  # A C G T
    return "".join(rng.choice(list("ACGT"), size=n, p=p))


# --------------------------------------------------------------------------
# 1. toy reference
# --------------------------------------------------------------------------
def build_reference(outdir, seed):
    rng = np.random.default_rng(seed)
    os.makedirs(outdir, exist_ok=True)
    contigs = {"chrT1": [], "chrT2": []}      # name -> list of sequence chunks
    gtf = []
    genes = {}                                 # name -> dict(strand, exons[(contig,start,end)], seq)
    pos = {c: 0 for c in contigs}

    def add(contig, seq):
        contigs[contig].append(seq)
        pos[contig] += len(seq)

    for i, name in enumerate(PANEL):
        contig = "chrT1" if i % 2 == 0 else "chrT2"
        strand = "+" if rng.random() < 0.6 else "-"
        add(contig, rand_seq(rng, int(rng.integers(500, 1500))))   # intergenic
        n_ex = int(rng.integers(3, 6))
        exons = []
        for k in range(n_ex):
            elen = int(rng.integers(160, 420))
            start = pos[contig] + 1                                  # GTF is 1-based
            seq = rand_seq(rng, elen)
            add(contig, seq)
            exons.append((start, start + elen - 1, seq))
            if k < n_ex - 1:
                add(contig, rand_seq(rng, int(rng.integers(300, 900))))  # intron
        g_start, g_end = exons[0][0], exons[-1][1]
        gid = f"TOYG{i + 1:011d}"
        tid = f"TOYT{i + 1:011d}"
        attr_g = (f'gene_id "{gid}"; gene_version "1"; gene_name "{name}"; '
                  f'gene_source "toy"; gene_biotype "protein_coding";')
        attr_t = (f'gene_id "{gid}"; gene_version "1"; transcript_id "{tid}"; '
                  f'transcript_version "1"; gene_name "{name}"; gene_source "toy"; '
                  f'gene_biotype "protein_coding"; transcript_name "{name}-201"; '
                  f'transcript_source "toy"; transcript_biotype "protein_coding";')
        gtf.append(f"{contig}\ttoy\tgene\t{g_start}\t{g_end}\t.\t{strand}\t.\t{attr_g}")
        gtf.append(f"{contig}\ttoy\ttranscript\t{g_start}\t{g_end}\t.\t{strand}\t.\t{attr_t}")
        order = range(len(exons)) if strand == "+" else range(len(exons) - 1, -1, -1)
        for rank, k in enumerate(order, 1):
            s, e, _ = exons[k]
            gtf.append(f"{contig}\ttoy\texon\t{s}\t{e}\t.\t{strand}\t.\t{attr_t} "
                       f'exon_number "{rank}"; exon_id "TOYE{i + 1:05d}{rank:02d}";')
        spliced = "".join(x[2] for x in exons)                      # + strand, mature
        mature = spliced if strand == "+" else revcomp(spliced)     # 5'->3' of the mRNA
        genes[name] = dict(gid=gid, mature=mature)
    for c in contigs:
        add(c, rand_seq(rng, 800))                                  # trailing spacer

    with open(os.path.join(outdir, "toy_genome.fa"), "w") as f:
        for c, chunks in contigs.items():
            seq = "".join(chunks)
            f.write(f">{c}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j:j + 60] + "\n")
    with open(os.path.join(outdir, "toy_genes.gtf"), "w") as f:
        f.write("\n".join(gtf) + "\n")
    tot = sum(len("".join(v)) for v in contigs.values())
    print(f"[reference] {len(PANEL)} genes, genome = {tot / 1e3:.0f} kb -> {outdir}")
    return genes


# --------------------------------------------------------------------------
# 2. toy FASTQ for one sample
# --------------------------------------------------------------------------
def load_real(sample_dir):
    """Return panel x spots count matrix (dense), in-tissue barcodes, off-tissue barcodes."""
    with h5py.File(os.path.join(sample_dir, "filtered_feature_bc_matrix.h5"), "r") as f:
        m = f["matrix"]
        names = [x.decode() for x in m["features/name"][:]]
        bcs = [x.decode().split("-")[0] for x in m["barcodes"][:]]   # drop "-1" GEM suffix
        data, ind, ptr = m["data"][:], m["indices"][:], m["indptr"][:]
        n_genes, n_spots = m["shape"][:]
    first = {}
    for i, n in enumerate(names):
        first.setdefault(n, i)
    missing = [g for g in PANEL if g not in first]
    if missing:
        sys.exit(f"panel genes not found in {sample_dir}: {missing}")
    rows = {first[g]: j for j, g in enumerate(PANEL)}
    mat = np.zeros((len(PANEL), n_spots), dtype=np.int64)     # CSC: columns = spots
    for s in range(n_spots):
        sl = slice(ptr[s], ptr[s + 1])
        for r, v in zip(ind[sl], data[sl]):
            j = rows.get(r)
            if j is not None:
                mat[j, s] = v
    pos = pd.read_csv(os.path.join(sample_dir, "spatial", "tissue_positions_list.csv"),
                      header=None, names=["bc", "in_tissue", "row", "col", "y", "x"])
    off = [b.split("-")[0] for b in pos.loc[pos.in_tissue == 0, "bc"]]
    return mat, bcs, off


def fq_writer(path):
    # mtime=0 and empty filename -> byte-reproducible gzip
    raw = open(path, "wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    return io.TextIOWrapper(gz, write_through=True), gz, raw


def make_sample_fastq(sample, sample_dir, genes, outdir, target_umis, seed,
                      background_frac=0.02, junk_frac=0.03, bc_error_frac=0.02):
    rng = np.random.default_rng(seed)
    pyr = random.Random(seed)
    mat, bcs, off_bcs = load_real(sample_dir)
    total = int(mat.sum())
    p = min(1.0, target_umis / total)
    thin = rng.binomial(mat, p)                                  # keep spatial pattern, fewer reads
    n_mol = int(thin.sum())

    # molecules: (barcode, gene, umi)
    gi, si = np.nonzero(thin)
    mols = []
    for g, s in zip(gi, si):
        mols.extend([(bcs[s], PANEL[g])] * int(thin[g, s]))
    # ambient background on off-tissue barcodes, gene frequencies of the sample
    n_bg = int(background_frac * n_mol)
    gfreq = thin.sum(axis=1) / thin.sum()
    for g in rng.choice(len(PANEL), size=n_bg, p=gfreq):
        mols.append((off_bcs[int(rng.integers(len(off_bcs)))], PANEL[g]))

    reads = []                                                    # (bc, umi, r2, kind)
    for bc, gname in mols:
        umi = rand_seq(rng, UMI_LEN, gc=0.5)
        mat_seq = genes[gname]["mature"]
        L = len(mat_seq)
        off = int(rng.integers(0, min(300, L - READ2_LEN) + 1))   # 3'-biased start
        r2 = mat_seq[L - READ2_LEN - off: L - off]
        n_dup = int(rng.choice([1, 2, 3, 4], p=[0.60, 0.25, 0.10, 0.05]))
        reads.extend([(bc, umi, r2, "gene")] * n_dup)
    n_junk = int(junk_frac * len(reads))                          # unmappable junk
    for _ in range(n_junk):
        reads.append((bcs[int(rng.integers(len(bcs)))], rand_seq(rng, UMI_LEN, 0.5),
                      rand_seq(rng, READ2_LEN), "junk"))
    pyr.shuffle(reads)

    os.makedirs(outdir, exist_ok=True)
    f1, g1, raw1 = fq_writer(os.path.join(outdir, f"{sample}_S1_L001_R1_001.fastq.gz"))
    f2, g2, raw2 = fq_writer(os.path.join(outdir, f"{sample}_S1_L001_R2_001.fastq.gz"))
    q1, q2 = "F" * READ1_LEN, "F" * READ2_LEN
    n_err = 0
    for i, (bc, umi, r2, _) in enumerate(reads):
        if rng.random() < bc_error_frac:                          # 1-mismatch barcode error
            k = int(rng.integers(BC_LEN))
            bc = bc[:k] + pyr.choice([b for b in "ACGT" if b != bc[k]]) + bc[k + 1:]
            n_err += 1
        name = f"@TOYSEQ:1:HTOY00DSX:1:{1101 + i // 20000}:{1000 + i % 20000}:{1000 + (i * 7) % 30000}"
        f1.write(f"{name} 1:N:0:ATCACGAT\n{bc}{umi}\n+\n{q1}\n")
        f2.write(f"{name} 2:N:0:ATCACGAT\n{r2}\n+\n{q2}\n")
    for f, g, raw in ((f1, g1, raw1), (f2, g2, raw2)):
        f.flush(); g.close(); raw.close()
    print(f"[{sample:7s}] {len(bcs):5d} in-tissue spots | UMIs {n_mol:6d} (+{n_bg} background) | "
          f"read pairs {len(reads):6d} ({n_junk} junk, {n_err} barcode errors)")
    return thin


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True,
                    help="GSE241124_RAW folder (contains Skin/, Wound1/, Wound7/, Wound30/)")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "toy_data"))
    ap.add_argument("--target-umis", type=int, default=25000,
                    help="approx. UMIs per sample after thinning (default 25000 -> ~40k read pairs)")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    genes = build_reference(os.path.join(a.out_dir, "reference_src"), a.seed)
    for i, s in enumerate(SAMPLES):
        make_sample_fastq(s, os.path.join(a.data_dir, s), genes,
                          os.path.join(a.out_dir, "fastq", s), a.target_umis, a.seed + 1 + i)
    print("done.")


if __name__ == "__main__":
    main()
