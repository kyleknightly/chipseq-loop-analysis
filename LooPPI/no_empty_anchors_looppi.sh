#!/usr/bin/env bash
set -euo pipefail

# This version filters OUT loops where NEITHER anchor has any proteins.
# Effected outputs (downstream computed from filtered loops):
#   - annotated_loops.bedpe         (replaced by filtered version)
#   - annotated_loop_ends.bed       (recomputed from filtered loops)
#   - loop_end_proportions.tsv      (recomputed)
#   - end_type_proportions.tsv      (recomputed)
#   - loop_type_proportions.tsv     (recomputed)
# If you also want cis/trans contact counts restricted, point those steps
# at the filtered loops as well (see notes near the bottom).

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
# Loop End IDs┬───────Anchor Details
#             │           │
#    Loop Details  Anchor Proportions
#    │    │    │
#    │    │   Loop End Proportions
#    │   End Type Proportions
#   Loop Type Proportions (using FILTERED loops)

LOOP_FILE="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/beds/hepg2-loops.bedpe"
CCRE_FILE="/mnt/altnas/work/Kyle.Knightly/chromatin-network/hepg2/beds/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.bed"
CHIP_DIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/inter-filtered"
OUT_DIR="./no-empty-out"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ====================================
printf "Creating anchor ID file \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
{
  k1=$1"_"$2"_"$3
  if(!(k1 in id)){ id[k1]=sprintf("A%06d", ++n); print $1,$2,$3,id[k1] }
  k2=$4"_"$5"_"$6
  if(!(k2 in id)){ id[k2]=sprintf("A%06d", ++n); print $4,$5,$6,id[k2] }
}' "$LOOP_FILE" > anchor_ids.bed

# ====================================
printf "Creating loop ID file \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{ map[$1"_"$2"_"$3]=$4; next }
{ print $1,$2,$3,map[$1"_"$2"_"$3],$4,$5,$6,map[$4"_"$5"_"$6] }' \
  anchor_ids.bed "$LOOP_FILE" > loop_ids.bedpe

# ====================================
printf "Finding anchor types \n" >&2
# ====================================
bedtools intersect -a anchor_ids.bed -b "$CCRE_FILE" -wa -wb -loj \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $14}' \
| sort -k4,4 -u > anchor_types.bed

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
  # Build per-id TF list (pipe-separated, unique)
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
printf "Calculating anchor proportions (unfiltered, for reference)\n" >&2
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
  for (tf in cnt) printf "%s\t%d\t%.6f\n", tf, cnt[tf], cnt[tf]/total
}' id_tfs.bed | sort -k1,1 > anchor_proportions.tsv

# ====================================
printf "Creating loop details file (with types & TFs per end) \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
FNR==NR { ty[$4]=$5; tf[$4]=$6; next }
{
  id1=$4; id2=$8
  print $1,$2,$3,id1,ty[id1],tf[id1],$5,$6,$7,id2,ty[id2],tf[id2]
}' annotated_anchors.bed loop_ids.bedpe > annotated_loops.bedpe

# ====================================
printf "FILTER: Drop loops with NO proteins on BOTH ends \n" >&2
# ====================================
# Keep a loop if (left has TFs) OR (right has TFs). Blank/NA/.'s are treated as empty.
awk 'BEGIN{FS=OFS="\t"}
{
  t1=$6; t2=$12
  if (t1=="." || t1=="NA") t1=""
  if (t2=="." || t2=="NA") t2=""
  if (t1!="" || t2!="") print $0
}' annotated_loops.bedpe > annotated_loops.filtered.bedpe

# Replace the original for downstream steps that you wanted impacted
mv -f annotated_loops.filtered.bedpe annotated_loops.bedpe

# ====================================
printf "Calculating loop end proportions (FILTERED) \n" >&2
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
  for (tf in cnt) printf "%s\t%d\t%.6f\n", tf, cnt[tf], (total? cnt[tf]/total : 0)
}' annotated_loop_ends.bed | sort -k1,1 > loop_end_proportions.tsv

# ====================================
printf "Calculating end type proportions (FILTERED) \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
{
  total++; if ($5!="") cnt[$5]++
}
END{
  for (t in cnt) printf "%s\t%d\t%.6f\n", t, cnt[t], (total? cnt[t]/total : 0)
}' annotated_loop_ends.bed | sort -k1,1 > end_type_proportions.tsv

# ====================================
printf "Calculating loop type proportions (FILTERED) \n" >&2
# ====================================
awk 'BEGIN{FS=OFS="\t"}
{
  total++                      # count every (filtered) loop row
  t1=$5; t2=$11
  if (t2<t1) { tmp=t1; t1=t2; t2=tmp }  # unordered pair
  pair[t1 SUBSEP t2]++
}
END{
  for (k in pair) {
    split(k,x,SUBSEP)
    printf "%s\t%s\t%d\t%.6f\n", x[1], x[2], pair[k], pair[k]/total
  }
}' annotated_loops.bedpe | sort -k1,1 -k2,2 > loop_type_proportions.tsv

# ====================================
# NOTE: cis/trans contact counting sections below still operate on the
# *unfiltered* anchor TF lists by default (id_tfs.bed or annotated_loops.bedpe).
# If you want those restricted to the filtered loop set as well, simply
# change their inputs from the original files to the filtered
# `annotated_loops.bedpe` produced above.
# ====================================

# (Optional) Recompute trans contacts on FILTERED loops
{
  awk 'BEGIN{FS=OFS="\t"}
  {
    split($6, A, /\|/); split($12, B, /\|/);
    delete SA; delete SB
    for (i in A) if (A[i]!="") SA[A[i]]=1
    for (j in B) if (B[j]!="") SB[B[j]]=1
    for (p in SA) for (q in SB) {
      if (p=="" || q=="") continue
      x=p; y=q; if (y<x) { tmp=x; x=y; y=tmp }
      pair[x SUBSEP y]++
    }
  }
  END{
    for (k in pair) { split(k, t, SUBSEP); print t[1], t[2], pair[k] }
  }' annotated_loops.bedpe | sort -k1,1 -k2,2
} > trans-contacts.tsv

printf "Done. Filtered loops written to annotated_loops.bedpe and downstream tables updated.\n" >&2
