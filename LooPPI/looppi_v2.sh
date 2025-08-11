#!/usr/bin/env bash

# ====================================
#              SCHEMATIC
# ====================================
# Loop list
#     │
#  End IDs
#     │
# Anchor IDs──┬──cCREs
#     │       │
#     │   Anchor Types────┬─ChIP-seq
#     │                   │
#  Loop IDs─┬───────Anchor Details
#           │             │
#    Loop Details  Anchor Proportions
#    │    │    │
#    │    │   Loop End Proportions
#    │   End Type Proportions
#   Loop Type Proportions

set -euo pipefail

LOOP_FILE="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe"
CCRE_FILE="/mnt/altnas/work/Kyle.Knightly/contact-network/hepg2/beds/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.bed"
CHIP_DIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/inter-filtered"
OUT_DIR="./out"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ====================================
printf "Creating anchor ID file \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="	"}
{
  k1=$1"_"$2"_"$3
  if(!(k1 in id)){ id[k1]=sprintf("A%06d", ++n); print $1,$2,$3,id[k1] }
  k2=$4"_"$5"_"$6
  if(!(k2 in id)){ id[k2]=sprintf("A%06d", ++n); print $4,$5,$6,id[k2] }
}' "$LOOP_FILE" > anchor_ids.bed

# ====================================
printf "Creating loop ID file \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="	"}
FNR==NR{ map[$1"_"$2"_"$3]=$4; next }
{ print $1,$2,$3,map[$1"_"$2"_"$3],$4,$5,$6,map[$4"_"$5"_"$6] }' \
  anchor_ids.bed "$LOOP_FILE" > loop_ids.bedpe

# ====================================
printf "Finding anchor types \n" >&2
# ====================================
bedtools intersect -a anchor_ids.bed -b "$CCRE_FILE" -wa -wb -loj \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $14}' | sort -k4,4 -u > anchor_types.bed

# ====================================
printf "Finding anchor proteins \n" >&2
# ====================================
> "hits.tmp"
for chip in "$CHIP_DIR"/*.bed; do
    tf_name=$(basename "$chip" | cut -d'-' -f1)
    echo "Processing $tf_name..." >&2

    bedtools intersect -a anchor_types.bed -b "$chip" -wa -u \
    | awk -v TF="$tf_name" 'BEGIN{OFS="\t"} {print $4, TF}' >> hits.tmp
done

sort -k1,1 -k2,2 -u hits.tmp > unique_tf_hits.tsv

awk 'BEGIN{FS=OFS="\t"}
FNR==NR {
  # Build per-id TF list
  if ($2!="") t[$1] = (t[$1] ? t[$1] "|" $2 : $2)
  next
}
{ print $4, t[$4] }' unique_tf_hits.tsv anchor_ids.bed > id_tfs.bed

awk '
BEGIN{FS=OFS="\t"}
FNR==NR {hits[$1]=$2; next}
{ print $1, $2, $3, $4, $5, hits[$4] }
' id_tfs.bed anchor_types.bed > annotated_anchors.bed

rm -f hits.tmp

# ====================================
printf "Calculating anchor proportions \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
{
  total++
  m = split($2, a, /\|/)
  delete seen
  for (i=1; i<=m; i++) {
    tf = a[i]
    if (tf != "" && !seen[tf]++) cnt[tf]++
  }
}
END{
  print "protein","count","proportion"
  for (tf in cnt) printf "%s\t%d\t%.6f\n", tf, cnt[tf], cnt[tf]/total
}' id_tfs.bed | sort -k1,1 > anchor_proportions.tsv

# ====================================
printf "Creating loop details file \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
FNR==NR { ty[$4]=$5; tf[$4]=$6; next }
{
  id1=$4; id2=$8
  print $1,$2,$3,id1,ty[id1],tf[id1],$5,$6,$7,id2,ty[id2],tf[id2]
}' annotated_anchors.bed loop_ids.bedpe > annotated_loops.bedpe

# ====================================
printf "Calculating loop end proportions \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}{
  print $1,$2,$3,$4,$5,$6;              # left end
  print $7,$8,$9,$10,$11,$12            # right end
}' annotated_loops.bedpe > annotated_loop_ends.bed
awk 'BEGIN{FS=OFS="\t"}
{
  total++
  m = split($6, a, /\|/)
  delete seen
  for (i=1; i<=m; i++) {
    tf = a[i]
    if (tf != "" && !seen[tf]++) cnt[tf]++
  }
}
END{
  print "protein","count","proportion"
  for (tf in cnt) printf "%s\t%d\t%.6f\n", tf, cnt[tf], cnt[tf]/total
}' annotated_loop_ends.bed | sort -k1,1 > loop_end_proportions.tsv
# ====================================
printf "Calculating end type proportions \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
{
  total++; cnt[$5]++
}
END{
  print "type","count","proportion"
  for (t in cnt) printf "%s\t%d\t%.6f\n", t, cnt[t], cnt[t]/total
}' annotated_loop_ends.bed | sort -k1,1 > end_type_proportions.tsv