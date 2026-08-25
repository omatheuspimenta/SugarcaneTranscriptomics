#!/bin/bash
set -euo pipefail

# ---- CONFIG ----
BAM_DIR="/home/diego/rnaseq/resultados_rnaseq/star_salmon"
STRINGTIE_DIR="/home/diego/rnaseq/resultados_rnaseq/star_salmon/stringtie"
REF_GTF="/dados04/jorge/diego_transcriptome/data/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf"
OUT_DIR="/dados04/jorge/diego_transcriptome/data/stringtie_denovo"
THREADS=200

mkdir -p "$OUT_DIR"
MERGELIST="$OUT_DIR/mergelist_denovo.txt"
> "$MERGELIST"   # empty/create the file

# ---- LOOP OVER SAMPLES ----
for gtf in "$STRINGTIE_DIR"/*.transcripts.gtf; do
    sample=$(basename "$gtf" .transcripts.gtf)
    bam="$BAM_DIR/${sample}.markdup.sorted.bam"

    if [[ ! -f "$bam" ]]; then
        echo "WARNING: BAM not found for sample $sample at $bam -- skipping" >&2
        continue
    fi

    out_gtf="$OUT_DIR/${sample}.denovo.transcripts.gtf"

    echo ">>> Running StringTie (de novo) on $sample"
    # If the read are stranded, you can add the -rf or -fr option to the command below. Check the StringTie manual for more details.
    stringtie "$bam" \
        -G "$REF_GTF" \
        -o "$out_gtf" \
        -p "$THREADS"

    echo "$out_gtf" >> "$MERGELIST"
done

echo "Done. Mergelist written to $MERGELIST"