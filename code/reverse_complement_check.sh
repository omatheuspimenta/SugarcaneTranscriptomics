#!/usr/bin/env bash
#
# reverse_complement_check.sh
#
# Bidirectional consistency check for lncRNA candidates that StringTie could not
# assign a strand to ('.' in column 7 of the GTF).
#
# Rationale
# ---------
# The libraries are unstranded, so monoexonic transcripts have no splice junction
# from which StringTie can infer orientation. `gffread -w` extracts these in the
# '+' orientation regardless, meaning that roughly half of them were submitted to
# the coding-potential predictors as the reverse complement of the real
# transcript. A monoexonic mRNA read backwards tends to score as non-coding, so
# this subset is enriched for false positives.
#
# This script classifies the same sequences a second time in the reverse
# orientation. Only transcripts predicted as lncRNA in BOTH orientations are
# retained.
#
# The consensus for the reverse orientation therefore uses CPC2, RNAsamba, RNAplonc and
# LncADeep2 (4 methods), against 5 for the forward orientation. Majority voting
# thresholds in consensus.R must be adjusted accordingly.
#
# Requires: seqkit, and the three predictors installed as described in README.md
#
set -euo pipefail

# Paths — adjust to your installation
#####################################
DATA=data
CPC2_BIN=/code/lncRNA/CPC2/CPC2_standalone-1.0.1/bin
RNASAMBA_DIR=/code/lncRNA/RNAsamba
LNCADEEP_DIR=/code/lncRNA/LncADeep2

REV_DIR="${DATA}/rev_predictions"
mkdir -p "${REV_DIR}"

#####################################
# List the transcripts with no strand assignment
#####################################
awk -F'\t' '$0!~/^#/ && $7=="."' "${DATA}/lncRNA_final_tagged.gtf" \
  | grep -o 'transcript_id "[^"]*"' \
  | sed 's/transcript_id "\(.*\)"/\1/' \
  | sort -u > "${DATA}/lncRNA_sem_fita.txt"

echo "Transcripts without strand: $(wc -l < "${DATA}/lncRNA_sem_fita.txt")"   # 1396

# confirm they are all monoexonic (expected: every transcript has 1 exon)
awk -F'\t' 'NR==FNR{ids[$1];next} $3=="exon"{
    match($9,/transcript_id "([^"]+)"/,a); if(a[1] in ids) n[a[1]]++
  } END{for(t in n) print n[t]}' \
  "${DATA}/lncRNA_sem_fita.txt" "${DATA}/lncRNA_final_tagged.gtf" \
  | sort | uniq -c

#####################################
# Extract those sequences and build their reverse complement
#####################################
seqkit grep -f "${DATA}/lncRNA_sem_fita.txt" "${DATA}/lncRNA_candidates.fa" \
  > "${DATA}/sem_fita.fa"

seqkit seq -r -p "${DATA}/sem_fita.fa" \
  | sed '/^>/s/$/_rev/' \
  > "${DATA}/sem_fita_rev.fa"

# sanity: both files must hold the same number of sequences as the ID list
grep -c '^>' "${DATA}/sem_fita.fa"
grep -c '^>' "${DATA}/sem_fita_rev.fa"

#####################################
# Re-run the predictors on the reverse-complemented sequences
#####################################

# --- CPC2 ---
( cd "${CPC2_BIN}" && \
  python CPC2.py \
    -i "$(realpath "${DATA}/sem_fita_rev.fa")" \
    -o "$(realpath "${REV_DIR}")/cpc2_rev" )

# --- RNAsamba ---
# uses the same pre-trained weights as the forward run
( cd "${RNASAMBA_DIR}" && \
  rnasamba classify \
    "$(realpath "${REV_DIR}")/rnasamba_rev.tsv" \
    "$(realpath "${DATA}/sem_fita_rev.fa")" \
    full_length_weights.hdf5 \
    partial_length_weights.hdf5 )

# --- LncADeep2 ---
( cd "${LNCADEEP_DIR}" && \
  python LncADeep2.py \
    -m identify \
    -i "$(realpath "${DATA}/sem_fita_rev.fa")" \
    -o "$(realpath "${REV_DIR}")/lncadeep2_rev.csv" \
    -d 'cuda:0' )

# --- RNAplonc ---
docker run --rm -v "$(realpath ${DATA})":/data rnaplonc \
    -i /data/sem_fita_rev.fa -o /data/rev_predictions/rnaplonc_rev.txt

#####################################
# Consensus over the reverse orientation
#####################################
# Run consensus.R over the three tables in ${REV_DIR}. It writes the IDs called
# lncRNA in the reverse orientation (still carrying the _rev suffix) to:
#   ${DATA}/lncRNA_rev_aprovados_raw.txt
Rscript code/consensus.R --mode reverse --indir "${REV_DIR}" \
  --out "${DATA}/lncRNA_rev_aprovados_raw.txt"

#####################################
# Intersect with the forward call and split approved / discarded
#####################################
sed 's/_rev$//' "${DATA}/lncRNA_rev_aprovados_raw.txt" | sort -u \
  > "${DATA}/lncRNA_semfita_aprovados.txt"

# every approved ID must belong to the unstranded set — must print 0
comm -23 "${DATA}/lncRNA_semfita_aprovados.txt" \
         <(sort -u "${DATA}/lncRNA_sem_fita.txt") | wc -l

comm -13 "${DATA}/lncRNA_semfita_aprovados.txt" \
         <(sort -u "${DATA}/lncRNA_sem_fita.txt") \
  > "${DATA}/lncRNA_semfita_reprovados.txt"

echo "Retained : $(wc -l < "${DATA}/lncRNA_semfita_aprovados.txt")"    # 670
echo "Discarded: $(wc -l < "${DATA}/lncRNA_semfita_reprovados.txt")"   # 726

#####################################
# Emit the retained transcripts with an arbitrary '+' strand, flagged
#####################################
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR { ok[$1]; next }
  $7=="." {
    match($9, /transcript_id "([^"]+)"/, a)
    if (a[1] in ok) { $7="+"; $9=$9" strand_arbitrary \"true\";"; print }
  }' "${DATA}/lncRNA_semfita_aprovados.txt" "${DATA}/lncRNA_final_tagged.gtf" \
  > "${DATA}/lncRNA_semfita_aprovados.gtf"

# 2 lines per transcript (one 'transcript', one 'exon') — expected 1340
wc -l "${DATA}/lncRNA_semfita_aprovados.gtf"
cut -f3 "${DATA}/lncRNA_semfita_aprovados.gtf" | sort | uniq -c
