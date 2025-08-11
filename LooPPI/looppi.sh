#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# loop_anchor_pipeline.sh  (multiplicity-weighted, single-intersect version)
# -----------------------------------------------------------------------------
# • Build UNIQUE anchors from BEDPE.
# • Intersect TF BEDs ONCE with the unique anchors.
# • Store anchor multiplicities (how many loop ends map to each unique anchor).
# • Compute TF proportions on unique anchors and multiplicity-weighted ends.
# • Build CIS matrix weighted by multiplicity; TRANS matrix per loop.
# -----------------------------------------------------------------------------
# Requirements: bash>=4, gawk, bedtools, coreutils, (optional) GNU parallel
# -----------------------------------------------------------------------------

set -euo pipefail

# ---- Defaults ----
DEFAULT_LOOPS="/mnt/altnas/work/Kyle.Knightly/looppi/random-hepg2/random-hepg2-loops.bedpe"
DEFAULT_PROT_DIR="/mnt/altnas/work/Kyle.Knightly/chipseq-analysis/hepg2/all-chipseq/new-merged-filtered"
DEFAULT_OUTDIR="./out"
DEFAULT_THREADS=1

# ---- Args ----
LOOPS_BEDPE=${1:-$DEFAULT_LOOPS}
PROT_DIR=${2:-$DEFAULT_PROT_DIR}
OUTDIR=${3:-$DEFAULT_OUTDIR}
THREADS=${4:-$DEFAULT_THREADS}

first_arg=${1-}
if [[ $first_arg == "-h" || $first_arg == "--help" ]]; then
  echo "Usage: $0 [loops.bedpe] [protein_bed_dir] [outdir] [threads]" >&2
  echo "Defaults:" >&2
  echo "  LOOPS_BEDPE   = $DEFAULT_LOOPS" >&2
  echo "  PROT_DIR      = $DEFAULT_PROT_DIR" >&2
  echo "  OUTDIR        = $DEFAULT_OUTDIR" >&2
  echo "  THREADS       = $DEFAULT_THREADS" >&2
  exit 0
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ---- Helpers ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found" >&2; exit 1; }; }
need bedtools; need gawk; need sort; need uniq

# -----------------------------------------------------------------------------
# [1/6] Explode BEDPE → NON-UNIQUE ends (2 per loop)
# -----------------------------------------------------------------------------
printf "[1/6] Exploding BEDPE into ends...\n" >&2
awk 'BEGIN{OFS="\t"}{
  id=sprintf("E%06d", NR);
  print $1,$2,$3,id"A";
  print $4,$5,$6,id"B";
}' "$LOOPS_BEDPE" > 00_ends.nonunique.bed

N_NONUNIQ=$(wc -l < 00_ends.nonunique.bed)

# -----------------------------------------------------------------------------
# [2/6] Collapse to UNIQUE anchors + anchor_id
# -----------------------------------------------------------------------------
printf "[2/6] Collapsing to unique anchors...\n" >&2
sort -k1,1 -k2,2n -k3,3n 00_ends.nonunique.bed \
| awk 'BEGIN{OFS="\t"}{k=$1"\t"$2"\t"$3; if(!(k in seen)){seen[k]=1; id=sprintf("A%06d", ++n); print $1,$2,$3,id}}' \
> 01_anchors.unique.bed

N_UNIQ=$(wc -l < 01_anchors.unique.bed)

# coord_key -> anchor_id map
awk 'BEGIN{OFS="\t"}{print $1"_"$2"_"$3, $4}' 01_anchors.unique.bed > anchors.map

# -----------------------------------------------------------------------------
# [3/6] Intersect UNIQUE anchors with protein BEDs (once)
# -----------------------------------------------------------------------------
printf "[3/6] Intersecting TF beds with UNIQUE anchors...\n" >&2

mapfile -t TF_BEDS < <(find "$PROT_DIR" -type f -name "*.bed" | sort)
if [[ ${#TF_BEDS[@]} -eq 0 ]]; then
  echo "ERROR: No *.bed files found in $PROT_DIR" >&2; exit 1
fi

TF_NAME_CMD='base=$(basename "$f"); base=${base%.bed}; echo ${base%%-*}'

> 02_anchor_tf_hits.unique.tsv

intersect_unique() {
  local f="$1"
  local tf
  tf=$(eval "$TF_NAME_CMD")
  bedtools intersect -a 01_anchors.unique.bed -b "$f" -wa -u \
    | awk -v TF="$tf" 'BEGIN{OFS="\t"}{print $4, TF}'
  echo "[3/6] done $tf" >&2
}

if (( THREADS > 1 )) && command -v parallel >/dev/null 2>&1; then
  export -f intersect_unique
  export TF_NAME_CMD OUTDIR
  parallel -j "$THREADS" intersect_unique ::: "${TF_BEDS[@]}" >> 02_anchor_tf_hits.unique.tsv
else
  for f in "${TF_BEDS[@]}"; do intersect_unique "$f" >> 02_anchor_tf_hits.unique.tsv; done
fi

sort -k1,1 -k2,2 02_anchor_tf_hits.unique.tsv -o 02_anchor_tf_hits.unique.tsv

# Anchor -> TF list
printf "[3a/6] Building per-anchor TF lists...\n" >&2
awk 'BEGIN{FS=OFS="\t"}{ tf_list[$1] = tf_list[$1] (tf_list[$1] ? "|" : "") $2 } END{ for(a in tf_list) print a, tf_list[a] }' \
  02_anchor_tf_hits.unique.tsv > anchor_tf_list.tsv

# -----------------------------------------------------------------------------
# [4/6] Anchor multiplicities (how many ends map to each anchor)
# -----------------------------------------------------------------------------
printf "[4/6] Computing anchor multiplicities...\n" >&2
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{coord2anchor[$1]=$2; next}
{
  key=$1"_"$2"_"$3
  aid=coord2anchor[key]
  if(aid!="") {mult[aid]++; total++}
}
END{
  for(a in mult) print a, mult[a] > "01a_anchor_multiplicity.tsv"
  print total > "TOTAL_END_MULT.txt"
}' anchors.map 00_ends.nonunique.bed

# -----------------------------------------------------------------------------
# [5/6] Protein proportions (unique + weighted)
# -----------------------------------------------------------------------------
printf "[5/6] Computing TF proportions (unique + weighted)...\n" >&2

# 5a unique
{
  echo -e "protein\tcount\tproportion"
  awk -F'\t' -v OFS='\t' -v N="$N_UNIQ" '
  { c[$2]++ }
  END{
    for (t in c) printf "%s\t%d\t%.6f\n", t, c[t], c[t]/N
  }' 02_anchor_tf_hits.unique.tsv | sort -k1,1
} > 04_tf_props.unique.tsv

# 5b weighted by multiplicity (non-unique ends surrogate)
{
  echo -e "protein\tcount\tproportion"
  awk -F'\t' -v OFS='\t' '
  FNR==NR { mult[$1]=$2; total+=mult[$1]; next }
  { c[$2] += (mult[$1] + 0) }
  END{
    for (t in c) printf "%s\t%d\t%.6f\n", t, c[t], c[t]/total
  }' 01a_anchor_multiplicity.tsv 02_anchor_tf_hits.unique.tsv | sort -k1,1
} > 05_tf_props.weightedEnds.tsv

# -----------------------------------------------------------------------------
# [6/6] CIS matrix (weighted) + TRANS matrix (per-loop)
# -----------------------------------------------------------------------------
printf "[6/6] Building CIS + TRANS matrices...\n" >&2

# CIS long (unweighted)
{
  echo -e "tf1\ttf2\tcount"
  awk 'BEGIN{FS=OFS="\t"}{
    n=split($2, arr, /\|/)
    for(i=1;i<=n;i++){
      for(j=i;j<=n;j++){
        tf1=arr[i]; tf2=arr[j]
        if(tf2<tf1){tmp=tf1; tf1=tf2; tf2=tmp}
        count[tf1,tf2]++
        seen[tf1]=1; seen[tf2]=1
      }
    }
  } END{
    for(k in count){split(k, idx, SUBSEP); print idx[1], idx[2], count[k]}
  }' anchor_tf_list.tsv | sort -k1,1 -k2,2
} > 06_cis_pairs.long.tsv

# CIS matrix
awk 'BEGIN{FS=OFS="\t"}
{pair[$1,$2]=$3; tf[$1]=1; tf[$2]=1}
END{
  i=0; for(t in tf){T[++i]=t}; asort(T)
  printf "TF"; for(j=1;j<=i;j++) printf "\t%s", T[j]; print ""
  for(r=1;r<=i;r++){
    a=T[r]; printf "%s", a
    for(c=1;c<=i;c++){
      b=T[c]; if(b<a){k1=b;k2=a}else{k1=a;k2=b}
      printf "\t%d", (pair[k1,k2]+0)
    }
    print ""
  }
}' 06_cis_pairs.long.tsv > 07_cis_matrix.tsv

# TRANS
awk 'BEGIN{FS=OFS="\t"}
FNR==NR{anchor_tf[$1]=$2; split($2, arr, /\|/); for(i in arr) all[arr[i]]=1; next}
ARGIND==2{coord2anchor[$1]=$2; next}
ARGIND==3{
  k1=$1"_"$2"_"$3; k2=$4"_"$5"_"$6
  a1=coord2anchor[k1]; a2=coord2anchor[k2]
  if(a1==""||a2=="") next
  t1=anchor_tf[a1]; t2=anchor_tf[a2]
  if(t1==""||t2=="") next
  n1=split(t1,A,/\|/); n2=split(t2,B,/\|/)
  for(i=1;i<=n1;i++) for(j=1;j<=n2;j++){
    tf1=A[i]; tf2=B[j]; if(tf2<tf1){tmp=tf1;tf1=tf2;tf2=tmp}
    pair[tf1,tf2]++
  }
}
END{
  print "tf1\ttf2\tcount" > "08_trans_pairs.long.tsv"
  for(k in pair){split(k,idx,SUBSEP); print idx[1], idx[2], pair[k] > "08_trans_pairs.long.tsv"}
  close("08_trans_pairs.long.tsv")
  i=0; for(t in all){T[++i]=t}; asort(T)
  printf "TF" > "09_trans_matrix.tsv"
  for(j=1;j<=i;j++) printf "\t%s", T[j] > "09_trans_matrix.tsv"
  printf "\n" > "09_trans_matrix.tsv"
  for(r=1;r<=i;r++){
    a=T[r]; printf "%s", a > "09_trans_matrix.tsv"
    for(c=1;c<=i;c++){
      b=T[c]; if(b<a){k1=b;k2=a}else{k1=a;k2=b}
      printf "\t%d", (pair[k1,k2]+0) > "09_trans_matrix.tsv"
    }
    printf "\n" > "09_trans_matrix.tsv"
  }
}' anchor_tf_list.tsv anchors.map "$LOOPS_BEDPE"

printf "Done. Outputs in: %s\n" "$OUTDIR" >&2
exit 0