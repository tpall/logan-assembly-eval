#!/bin/bash
# Pathracer L1 demo -- why graph-aware recovery matters for Logan, and the exact
# recipe from Shen et al., "Petabase-scale Papillomavirus Discovery" (bioRxiv
# 2026.04.21.719858; code: github.com/syueqiao/pvs_in_logan).
#
# It shows the PROBLEM live (papillomavirus L1 shattered across Logan contigs)
# and reconstructs the de Bruijn graph from the contig headers -- the substrate
# graph-aware recovery runs on. The RECOVERY step uses `logan-extract` +
# `pathracer`, which are NOT on bioconda: they come from the SPAdes `graph-extract`
# branch (github.com/ablab/spades/tree/graph-extract) + PathRacer in SPAdes v4.0.0.
# If those two tools are on PATH the demo runs the real recipe; otherwise it prints
# it. (A generic bioconda PathRacer will NOT work -- see demo/README.md.)
#
# Deps (baseline): hmmer prodigal seqkit zstd curl + awk.  Run: bash run_l1_demo.sh [ACC]
set -euo pipefail
ACC=${1:-SRR2565980}          # HPV-rich metagenome (~1e5 HPV reads)
K=31                          # Logan contigs: Minia3, de Bruijn k=31
HMM=${HMM:-PF00500.hmm}       # Pfam papillomavirus L1 (the paper uses a bespoke jrHMM)
THREADS=${THREADS:-8}
cd "$(dirname "$0")"
[ -s "$HMM" ] || { [ -s "$HMM.gz" ] && gunzip -kf "$HMM.gz"; }

echo "== 1. Logan contigs for $ACC =="
if [ ! -s "$ACC.contigs.fa" ]; then
  curl -sS -o "$ACC.contigs.fa.zst" "https://s3.amazonaws.com/logan-pub/c/$ACC/$ACC.contigs.fa.zst"
  zstd -dq -f "$ACC.contigs.fa.zst"; rm -f "$ACC.contigs.fa.zst"
fi
echo "   $(grep -c '^>' "$ACC.contigs.fa") contigs"

echo "== 2. the graph is in the headers: rebuild a GFA (illustrative) =="
awk -v ovl="${K}M" '
  /^>/ { name=substr($1,2); sub(/^.*_/,"",name); cov="0"; delete lk; nl=0
    for(i=2;i<=NF;i++){ if($i~/^ka:f:/){split($i,a,":");cov=a[3]} else if($i~/^L:/){lk[++nl]=$i} }
    getline seq; printf "S\t%s\t%s\tLN:i:%d\tKC:f:%s\n", name, seq, length(seq), cov
    for(j=1;j<=nl;j++){ split(lk[j],p,":"); printf "L\t%s\t%s\t%s\t%s\t%s\n", name, p[2], p[3], p[4], ovl }
  }' "$ACC.contigs.fa" > "$ACC.gfa"
echo "   $(grep -c '^S' "$ACC.gfa") nodes, $(grep -c '^L' "$ACC.gfa") edges"
echo "   NOTE: this GFA is a compacted-dBG (overlap k-1); generic PathRacer wants a"
echo "   SPAdes-convention graph (overlap = odd k). logan-extract does the real rebuild."

echo "== 3. THE PROBLEM: L1 on the isolated contigs (prodigal + hmmsearch) =="
[ -s "$ACC.faa" ] || prodigal -p meta -q -i "$ACC.contigs.fa" -a "$ACC.faa" -d /dev/null -o /dev/null 2>/dev/null
[ -s "$ACC.baseline.domtbl" ] || hmmsearch --noali --domtblout "$ACC.baseline.domtbl" -E 1e-5 "$HMM" "$ACC.faa" >/dev/null 2>&1 || true
awk '$1!~/^#/{printf "     L1 hit: %-24s aa=%-5s env=%s-%s\n",$1,$3,$20,$21}' "$ACC.baseline.domtbl" | sort -u | head -10
echo "   -> a full ~500 aa L1 is scattered across contigs as partial hits: that is the"
echo "      fragmentation graph-aware recovery is designed to stitch back together."

echo "== 4. RECOVERY (Shen et al.): logan-extract (default k) -> PathRacer =="
LE=${LOGAN_EXTRACT:-logan-extract}; PR=${PATHRACER:-pathracer}
if command -v "$LE" >/dev/null && command -v "$PR" >/dev/null; then
  # seed list = hit names in "diamond" format; logan-extract keys edges off split("_")[1]
  awk '$1!~/^#/{print $1}' "$ACC.baseline.domtbl" | sort -u > "${ACC}_list.txt"
  mkdir -p tmp
  # IMPORTANT: no -k. The default (k=29) rebuilds a *linked* subgraph; -k 31 yields a
  # link-less GFA that PathRacer cannot read ("Failed to determine k-mer length").
  "$LE" "$ACC.contigs.fa" "${ACC}_list.txt" "$ACC.subgraph.gfa" -d 1000 -t "$THREADS" -tmp-dir tmp --gfa
  echo "   subgraph: $(grep -c '^S' "$ACC.subgraph.gfa") segments, $(grep -c '^L' "$ACC.subgraph.gfa") links"
  "$PR" "$HMM" "$ACC.subgraph.gfa" --hmm --length 0.95 --top 5 --rescore -t "$THREADS" --output "pr_$ACC"
  echo "   recovered full-length L1 (top paths; Edges=A_B_C = stitched across graph edges A,B,C):"
  grep '^>' "pr_$ACC/${HMMNAME:-Late_protein_L1}.seqs.fa" 2>/dev/null | sed -E 's/\|Scaffold.*//' | head -5
else
  cat <<'RECIPE'
   [logan-extract / pathracer not on PATH -- build them once from the SPAdes
    graph-extract branch (verified 2026-07-27), then re-run this script:]

     git clone -b graph-extract https://github.com/ablab/spades && cd spades
     # logan-extract builds via the implicit-projects path (needs zstd):
     PREFIX="$PWD/inst" ./spades_compile.sh -j 8            # -> inst/bin/logan-extract
     # pathracer is OFF by default -- enable it explicitly at a FRESH configure:
     mkdir bp && cd bp
     cmake -DSPADES_ENABLE_PROJECTS="pathracer" -DCMAKE_INSTALL_PREFIX="$PWD/../inst" ../src
     make -j 8 pathracer                                    # -> ./projects/pathracer/pathracer

   Then: LOGAN_EXTRACT=/path/to/logan-extract PATHRACER=/path/to/pathracer bash run_l1_demo.sh
   (Paper uses a bespoke jrHMM over the 6 L1 sheets; PF00500 is the public stand-in here.)
RECIPE
fi
echo "DONE"