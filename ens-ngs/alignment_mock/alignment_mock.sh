#!/usr/bin/env bash
# alignment_mock.sh -- mock Space Ranger alignment on TOY FASTQs (ENS NGS practical)
# The reference is assumed to exist already (built once by the instructor, see below).
# Results are meaningless (toy data): run it, look at the logs / web_summary.html,
# then go back to the notebook and use the real, already aligned data.

# --- instructor, once (not for students) -------------------------------------
# The FASTA/GTF in toy_data/reference_src are only the INPUT of mkref. What
# --transcriptome needs is the folder that mkref creates (toy_data/toy_ref).
# cd toy_data
# spaceranger mkref --genome=toy_ref \
#           --fasta=reference_src/toy_genome.fa \
#           --genes=reference_src/toy_genes.gtf
# cd ..

for SAMPLE in Skin Wound1 Wound7 Wound30; do
    spaceranger count --id="${SAMPLE}" \
              --transcriptome=toy_data/toy_ref \
              --fastqs=toy_data/fastq/${SAMPLE} \
              --sample=${SAMPLE} \
              --image=toy_data/images/${SAMPLE}.jpg \
              --unknown-slide=visium-1 \
              --reorient-images=true \
              --create-bam=false \
              --localcores=8 \
              --localmem=32
done
